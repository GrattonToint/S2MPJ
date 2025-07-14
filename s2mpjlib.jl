#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                                S2MPJ library for Julia
#
#   Performs the runtime actions specific to S2MPJ, irrespective of the problem at hand.
#
#   Programming: S. Gratton and Ph. L. Toint (this version 14 VII 2025)
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

using SparseArrays

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Declare the PB and PBM data structures
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mutable struct PB
    name::String
    sifpbname::String
    n::Int64
    nob::Int64
    nle::Int64
    neq::Int64
    nge::Int64
    m::Int64
    lincons::Vector{Int64}
    pbclass::String
    x0::Vector{Float64}
    xlower::Vector{Float64}
    xupper::Vector{Float64}
    xscale::Vector{Float64}
    xtype::Vector{String}
    y0::Vector{Float64}
    clower::Vector{Float64}
    cupper::Vector{Float64}
    objlower::Float64
    objupper::Float64
    xnames::Vector{String}
    cnames::Vector{String}
    objderlvl::Int64
    conderlvl::Vector{Int64}

    PB(name::String) = new( name,      #name
                            "",        #sifpbname
                            0,         #n
                            0,         #nob
                            0,         #nle
                            0,         #neq
                            0,         #nge
                            0,         #m
                            Int64[],   #lincons
                            "",        #pbclass
                            Float64[], #x0
                            Float64[], #xlower
                            Float64[], #xupper
                            Float64[], #xscale
                            String[],  #xtype
                            Float64[], #y0
                            Float64[], #clower
                            Float64[], #cupper
                            -Inf,      #objlower
                             Inf,      #objupper
                            String[],  #xnames
                            String[],  #cnames
                             2,        #objderlvl
                             Int64[]   #conderlvl
                            )
end

mutable struct PBM
    name::String
    objgrps::Vector{Int64}
    congrps::Vector{Int64}
    A::SparseMatrixCSC{Float64, Int64}
    gconst::Vector{Float64}
    H::SparseMatrixCSC{Float64, Int64}
    enames::Vector{String}
    elftype::Vector{String}
    elvar::Vector{Vector{Int64}}
    elpar::Vector{Vector{Float64}}
    gscale::Vector{Float64}
    grnames::Vector{String}
    grftype::Vector{String}
    grelt::Vector{Vector{Int64}}
    grelw::Vector{Vector{Float64}}
    grpar::Vector{Vector{Float64}}
    efpar::Vector{Float64}
    gfpar::Vector{Float64}
    has_globs::Vector{Int64}
    objderlvl::Int64
    conderlvl::Vector{Int64}
    call::Function
    
    PBM(name::String) = new( name,                             #name
                             Int64[],                          #objgrps
                             Int64[],                          #congrps
                             spzeros(Float64,0,0),             #A
                             Float64[],                        #gconst
                             spzeros(Float64,0,0),             #H
                             String[],                         #enames
                             String[],                         #elftype
                             Vector{Vector{Int64}}(),          #elvar
                             Vector{Vector{Float64}}(),        #elpar
                             Float64[],                        #gscale
                             String[],                         #grnames
                             String[],                         #grftype
                             Vector{Vector{Int64}}(),          #grelt
                             Vector{Vector{Float64}}(),        #grelw
                             Vector{Vector{Float64}}(),        #grpar
                             Float64[],                        #efpar
                             Float64[],                        #gfpar
                             [0, 0],                           #has_globs
                             2,                                #objderlvl
                             Int64[]                           #conderlvl
                             )
end


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#   Computes the effective index of name in List, or add name to List if not in there already.
#   Return the index of name in the dictionary and new = 1 if the dictionary has been enlarged
#   or 0 otherwise.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s2mpj_ii( name::String, dict::Dict{String,Int} )
    if haskey( dict, name )
        new = 0
        idx = dict[name]
    else
        new = 1
        idx = length(dict) + 1
        dict[name] = idx
    end
    return idx, dict, new
end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#   Get the index of a nonlinear variable.  This implies adding it to the variables' dictionary ix_
#   if it is a new one, with bounds, start point and types defined by their default settings (it is
#   too late to define problem-specific values).
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s2mpj_nlx( name::String, dict::Dict{String,Int}, pb::PB, getxnames::Int,
                  xlowdef::Union{Float64,Nothing}, xuppdef::Union{Float64,Nothing}, x0def::Union{Float64,Nothing} )

    iv, dict, newvar = s2mpj_ii( name, dict )

    if newvar == 1
        pb.n += 1
        if !isnothing( getxnames )
            arrset( pb.xnames, iv, name )
        end
        if hasproperty( pb, :xlower )
            arrset( pb.xlower, iv, isnothing( xlowdef ) ? -Inf : xlowdef )
        end
        if hasproperty( pb, :xupper )
            arrset( pb.xupper, iv, isnothing( xuppdef ) ? Inf : xuppdef )
        end
        try
            arrset( pb.xtype, iv, "r" )
        catch
            # Ignore any errors in setting xtype
        end
        arrset( pb.x0, iv, isnothing(x0def) ? 0.0 : x0def )
    end

    return iv, dict, pb
    
end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#   Insert a value in a list, possibly enlarging the list to do so.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function arrset( thelist::Vector{T}, iv::Int, value::T ) where T

    if iv > length( thelist )
       if isa( value, String )
           while length(thelist) < iv
               push!( thelist, "" )  # Add elements to reach iv2
           end
        elseif isa( value, Int64 )
           append!(thelist, zeros(Int,iv-length(thelist) ) )
        elseif isa( value, Float64 )
           append!( thelist, zeros(Float64,iv-length(thelist) ) )
        else
           resize!( thelist, iv )
        end
    end
    thelist[iv] = value
    
end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#   Insert a value in a list of lists, possibly enlarging the list(s) to do so.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function loaset( list::Vector{Vector{T}}, iv1::Int, iv2::Int, value::T ) where T
    while length(list) < iv1
        push!(list, Vector{T}())   # Add empty vectors to reach iv1
    end
    arrset( list[iv1], iv2, value )
end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#   Perform the main high-level action requested from a S2MPJ problem file, using information computed
#   during the 'setup' call and passed to this function in the pbm struct.  The nonlinear functions
#   are evaluated at x = args[2]. See the documentation of S2MPJ for details.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s2mpj_eval( action::String, pbm::PBM, args::Union{Vector{Int},Vector{Float64}}... )

    #  Check plausibility of the request.

    if action in  ["fx","fgx","fgHx","fHxv"]
       if isempty( pbm.objgrps ) && isempty( pbm.H )
          println(" ")
          println("ERROR: problem "*pbm.name*".jl has no objective function!")
          println("       Please refer to the problem classification for checking a problem's type.")
          println(" ")
          return Nothing
       end
    elseif action in  ["cx","cJx","cJHx","cIx","cIJx","cIJHx","cIJxv","cJxv","cJtxv","cIJtxv"]
       if isempty( pbm.congrps )
          println(" ")
          println("ERROR: problem "*pbm.name*".jl has no constraint!")
          println("       Please refer to the problem classification for checking a problem's type.")
          println(" ")
          return Nothing
       end
    end

    #  Get the problem data and possible global parameters' values.
    
    if pbm.has_globs[1] >0 
        pbm = pbm.call( "e_globs", pbm )
    end
    if pbm.has_globs[2] >0 
        pbm = pbm.call( "g_globs", pbm )
    end
    
    #  Evaluations of the objective function and constraints

    print(pbm.objgrps)
    
    if action == "fx" 
        return evalgrsum( true,  pbm.objgrps, args[1], pbm, 1 )
    elseif action == "fgx" 
        return evalgrsum( true,  pbm.objgrps, args[1], pbm, 2 )
    elseif action == "fgHx"
        return evalgrsum( true,  pbm.objgrps, args[1], pbm, 3 )
    elseif action == "cx"
        return evalgrsum( false, pbm.congrps, args[1], pbm, 1 )
    elseif action == "cJx"
        return evalgrsum( false, pbm.congrps, args[1], pbm, 2 )
    elseif action == "cJHx"
        return evalgrsum( false, pbm.congrps, args[1], pbm, 3 )
    elseif action == "cIx"
        return evalgrsum( false, pbm.congrps[args[1]], args[1], pbm, 1 )
    elseif action == "cIJx"
        return evalgrsum( false, pbm.congrps[args[1]], args[1], pbm, 2 )
    elseif action == "cIJHx"
        return evalgrsum( false, pbm.congrps[args[1]], args[1], pbm, 3 )
    elseif action == "fHxv"
        return evalHJv( "Hv", pbm.objgrps, args[1], args[2], Float64[], pbm ); 
    elseif action == "cJxv"
        return evalHJv( "Jv", pbm.congrps, args[1], args[2], Float64[], pbm );
    elseif action == "cJtxv"
        return evalHJv( "Jtv", pbm.congrps, args[1], args[2], Float64[], pbm );
    elseif action == "cIJxv"
        return evalHJv( "Jv", pbm.congrps[args[3]], args[1], args[2], Float64[], pbm )
    elseif action == "cIJtxv"
        return evalHJv( "Jtv", pbm.congrps[args[3]], args[1], args[2], Float64[], pbm )
   
    # For the Lagrangian
    
    elseif action == "Lxy"
        return evalLx( pbm.objgrps, pbm.congrps, args[1], args[2], pbm, 1 );
    elseif action == "Lgxy"
        return evalLx( pbm.objgrps, pbm.congrps, args[1], args[2], pbm, 2 );
    elseif action == "LgHxy"
        return evalLx( pbm.objgrps, pbm.congrps, args[1], args[2], pbm, 3 );
    elseif action == "LIxy"
        return evalLx( pbm.objgrps, pbm.congrps[args[3]], args[1], args[2], pbm, 1 )
    elseif action == "LIgxy"
        return evalLx( pbm.objgrps, pbm.congrps[args[3]], args[1], args[2], pbm, 2 )
    elseif action == "LIgHxy"
        return evalLx( pbm.objgrps, pbm.congrps[args[3]], args[1], args[2], pbm, 3 )
    elseif action == "LHxyv"
        return evalLHxyv( pbm.objgrps, pbm.congrps, args[1], args[2], args[3], pbm )
    elseif action == "LIHxyv"
        return evalLHxyv( pbm.objgrps, pbm.congrps[args[4]], args[1], args[2], args[3], pbm )

    # Error

    else
        println("ERROR: action "*action*" unavailable for problem "*pbm.name*".jl")
        return ntuple(i->undef,args[end])
    end
    
end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#   Evaluate the value of a sum of groups (and, if requested, that of of its gradient
#   and Hessian) at x, given the problem data available in the pbm struct.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function evalgrsum( isobj::Bool, glist::Vector{Int}, x::Vector{Float64}, pbm::PBM, nargout::Int )

   debug = false#D

#  Initializations

   n  = length( x )
   if isobj
      fx = 0
   else
      m    = length( glist )
      cx   = zeros( Float64, m )
      ic   = 0
      if !isempty(pbm.conderlvl)
         lder = length( pbm.conderlvl )
      end
   end
   if nargout > 1
      if isobj
         if !isempty(pbm.objderlvl)
             if pbm.objderlvl >= 1
                 gx = spzeros( n, 1 )
             else
                 gx = spzeros( n, 1 )
                 gx[1] = NaN
             end
          else
             gx = spzeros( n, 1 )
          end
      else
         if !isempty(pbm.conderlvl)
             if ( lder == 1 && pbm.conderlvl[1] >= 1) || lder > 1
                Jx = spzeros( m, n )
             else
                Jx = spzeros( m, n )
                Jx[1,1] = NaN
             end
         else
             Jx = spzeros( m, n )
         end         
      end
      if nargout > 2
         if isobj
            if !isempty(pbm.objderlvl)
               if pbm.objderlvl >= 2
                   Hx = spzeros( n, n )
               else
                   Hx = spzeros( n, n )
                   Hx[1,1] = NaN
               end
            else
               Hx = spzeros( n, n )
            end            
         else
            Hx = []
         end
      end
   end

   # Check if pbm.A is not empty

   if !isempty(pbm.A)
       sA1, sA2 = size(pbm.A)
       has_A = true
   else
       has_A = false
   end
   
   # Evaluate the quadratic term, if any.
   
   if isobj && !isempty(pbm.H)
       Htimesx = pbm.H * x
       if nargout == 1
           fx += 0.5 * x'* Htimesx
       elseif nargout == 2
           gx += Htimesx
           fx += 0.5 * x'* Htimesx
       elseif nargout == 3
           Htimesx = pbm.H * x # This line appears redundant and could be omitted.
           gx += Htimesx
           fx += 0.5 * x'* Htimesx
           Hx += pbm.H
       end
   end

   if debug
       if isobj
           println( "fx(quadratic) = ", fx )
       else
           println( "cx(quadratic) = ", cx )
       end
   end

   for iig in 1:length(glist)
       ig = glist[iig]
       #  Find the level of available derivatives for the group.

       if isobj
           if !isempty(pbm.objderlvl)
               derlvl = pbm.objderlvl
           else
               derlvl = 2
           end
       else
           if !isempty(pbm.conderlvl)
               if lder == 1
                   derlvl = pbm.conderlvl[ 1 ]
               else
                   derlvl = pbm.conderlvl[ findfirst( isequal(ig), pbm.congrps ) ]
               end
           else
               derlvl = 2
           end
       end
       nout = min( nargout, derlvl + 1 );

       # Find the group's scaling.
       
       gsc = 1.0 # Default value
       if isdefined(pbm, :gscale)
           if ig <= length(pbm.gscale) && abs(pbm.gscale[ig]) > 1.0e-15
               gsc = pbm.gscale[ig]
           end
       end

       # Evaluate the linear term, if any.

       if isdefined(pbm, :gconst) && ig <= length(pbm.gconst)
           fin = -pbm.gconst[ig]
       else
           fin = 0
       end
       if nargout == 1
           if has_A && ig <= sA1
               fin += pbm.A[ig, 1:sA2]' * x[1:sA2]
           end
       elseif nargout in [2, 3]
           gin = zeros( Float64, n )
           if has_A && ig <= sA1
               gin[1:sA2] = pbm.A[ig, 1:sA2]
               fin += gin[1:sA2]' * x[1:sA2]
           end
       end
   
       if nargout > 2
           Hin = spzeros(n, n)
       end
   
      if debug
          println("ig = $ig  fin(linear) = $fin")
      end
   
       # Loop on the group's elements.

       if isdefined(pbm, :grelt)  && ig <= length(pbm.grelt) && !isempty(pbm.grelt[ig])
           for iiel in 1:length(pbm.grelt[ig])
               iel    = pbm.grelt[ig][iiel]
               irange = pbm.elvar[iel]
               efname = pbm.elftype[iel]
               if isdefined( pbm, :grelw ) && ig <= length( pbm.grelw ) && !isempty( pbm.grelw[ig] )
                   has_weights = true
                   wiel        = pbm.grelw[ig][iiel]
               else
                   has_weights = false
               end
               if nout == 1
                   fiel = pbm.call( efname, x[irange], iel, 1, pbm )
                   fin += has_weights ? wiel * fiel : fiel
               elseif nout == 2
                   fiel, giel = pbm.call( efname, x[irange], iel, 2, pbm )
                   fin       += has_weights ? wiel * fiel : fiel
                   for ir in 1:length(irange)
                       ii       = irange[ir]
                       gin[ii] += has_weights ? wiel * giel[ir] : giel[ir]
                   end
               elseif nout == 3
                   fiel, giel, Hiel = pbm.call( efname, x[irange], iel, 3, pbm )
#                   println( "irange = ", irange, "  x[irange] = ", x[irange],
#                            "  efname = ", efname, "  iel = ", iel, "  fiel = ", fiel, " wiel = ", wiel )#D
#                   println(" Hiel = ", Hiel ) #D
                   fin += has_weights ? wiel * fiel : fiel
                   for ir in 1:length(irange)
                       ii       = irange[ir]
                       gin[ii] += has_weights ? wiel * giel[ir] : giel[ir]
                       for jr in 1:length(irange)
                          jj          = irange[jr]
                          Hin[ii,jj] += has_weights ? wiel * Hiel[ir,jr] : Hiel[ir,jr]
                       end
                   end
               end
           end
       end

       if debug
           if isobj
               println( "ig = ", ig, "  fin(nonlinear) = ", fin, "  fx = ", fx  )
           else
               println( "ig = ", ig, "  fin(nonlinear) = ", fin, "  cx = ", cx  )
           end
       end
       
       #  Evaluate the group function.
       #  1) the non-TRIVIAL case

       if isdefined(pbm, :grftype) && ig <= length(pbm.grftype) &&
          isassigned(pbm.grftype,ig) && pbm.grftype[ig] != "TRIVIAL" && pbm.grftype[ig] != ""
           egname = pbm.grftype[ig]  # the group's ftype
    
           if isobj
           
               if nargout == 1
                   fa = pbm.call( egname, fin, ig, 1, pbm )
                   fx += fa / gsc
               elseif nargout == 2
                   fa, grada = pbm.call( egname, fin, ig, 2, pbm )
                   fx += fa / gsc
                   if derlvl >= 1
                       gx += grada[1] * gin / gsc
                   else
                       gx = spzeros( n, 1 )
                       gx[1] = NaN
                   end
               elseif nargout == 3
#                   println( "fin = ", fin, " egname = ", egname, " gsc = ", gsc )#D
                   fa, grada, Hessa = pbm.call( egname, fin, ig, 3, pbm )
                   fx += fa / gsc
                   if derlvl >= 1
                       gx += grada[1] * gin / gsc
                   else
                       gx = spzeros( n, 1 )
                       gx[1] = NaN
                   end
                   if derlvl >= 2
                       sgin = sparse(gin)
                       Hx += (Hessa[1] * sgin * sgin' + grada[1] * Hin) / gsc
                   else
                       Hx = spzeros( n, n )
                       Hx[1,1] = NaN
                   end
               end
           else
               ic += 1
               if nargout == 1
                   fa = pbm.call( egname, fin, ig, 1, pbm )
                   cx[ic] = fa / gsc
               elseif nargout == 2
                   fa, grada = pbm.call( egname, fin, ig, 2, pbm )
                   cx[ic] = fa / gsc
                   if derlvl >= 1
                      Jx[ic, 1:n] = (grada[1] * gin') / gsc
                   else
                      Jx[ib, 1:n] = NaN* ones( 1, n )
                   end
               elseif nargout == 3
                   fa, grada, Hessa = pbm.call( egname, fin, ig, 3, pbm )
                   cx[ic] = fa / gsc
                   sgin = sparse(gin)
                   if derlvl >= 1
                      Jx[ic, 1:n] = (grada[1] * gin') / gsc
                   else
                      Jx[ib, 1:n] = NaN* ones( 1, n )
                   end
                   if derlvl >= 2
                      push!(Hx, (Hessa[1] * sgin * sgin' + grada[1] * Hin) / gsc)
                   else
                      Hin = spzeros( n, n )
                      Hin[1,1] = NaN
                      push!(Hx, Hin )
                   end
               end
           end
       else
       
           # TRIVIAL case
           
           if isobj
#               println( " objective fin = ", fin, " egname = TRIVIAL", " gsc = ", gsc )#D
               fx += fin / gsc
               if nargout >= 2
                   if derlvl >= 1
                      gx += gin / gsc
                   else
                      gx = spzeros( n, 1 )
                      gx[1] = NaN
                   end
               end
               if nargout == 3
                   if derlvl >= 2
                      Hx += Hin / gsc
                   else
                      Hx = spzeros( n, n )
                      Hx[1,1] = NaN
                   end
               end
           else
#               println( " constraint fin = ", fin, " egname = TRIVIAL", " gsc = ", gsc )#D
               ic += 1
               cx[ic] = fin / gsc
               if nargout >= 2
                   if derlvl >= 1
                      Jx[ic, 1:n] = gin' / gsc
                   else
                      Jx[ic, 1:n] = NaN * ones( 1, n )
                   end
               end
               if nargout == 3
                   if derlvl >= 2
                       push!( Hx, Hin / gsc)
                   else
                       Hin = spzeros( n, n )
                       Hin[1,1] = NaN
                       push!( Hx, Hin )
                   end
               end
           end
       end

#       println( "ig = ", ig, " Hx(final) = ", Hx )#D
       if debug
            if isobj
                println( "ig = ", ig, "  fx(final) = ", fx )
            else
                println( "ig = ", ig, "  cx[ic](final) = ", cx[ic] )
            end
        end
    
   end
   
    # Initialize output based on nargout

    varargout = Tuple{}
    if nargout == 1
        varargout = isobj ? fx : cx
    elseif nargout > 1
        varargout = isobj ? (fx,) : (cx,)
        additional_output = isobj ? gx : Jx
        varargout = (varargout..., additional_output)
        if nargout > 2
            varargout = (varargout..., Hx) # Hx is appended regardless of isobj; 
        end
    end
    
    return varargout

end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#   Evaluate the product of the Hessian of a sum of groups at x times a user-supplied vector v,
#   given the problem data available in the pbm struct.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function evalHJv( mode::String, glist::Vector{Int}, x::Vector{Float64},
                  v::Vector{Float64}, y::Vector{Float64}, pbm::PBM )

    debug = false#D
   
    #  Initializations
    #  Avoid computing anything when the relevant derivatives are missing.

    n  = length( x );
    if mode == "Hv"
        if !isempty(pbm.objderlvl)
           if ( pbm.objderlvl < 2 )
              HJv    = zeros( Float64, n )
              HJv[1] = NaN
              return HJv
           else
              HJv    = zeros( Float64, n )
              derlvl = 2
           end
        else
           HJv    = zeros( Float64, n )
           derlvl = 2
        end        
    elseif mode == "HIv"
        m = length( glist )
        if !isempty(pbm.conderlvl)
            if ( any( x -> x<2, pbm.conderlvl ) )
                HJv    = zeros( Float64, n )
                HJv[1] = NaN
                return HJv
            else
                HJv    = zeros( Float64, n )
                ic     = 0
                derlvl = 2 
            end
        else
            HJv    = zeros( Float64, n )
            ic     = 0
            derlvl = 2 
        end
    else
        m = length( glist )
        if !isempty(pbm.conderlvl)
            if ( any( x -> x<1, pbm.conderlvl ) )
                if mode == "Jv"
                    HJv = zeros( Float64, m )
                else
                    HJv = zeros( Float64, n )
                end
                HJv[1] = NaN
                return HJv
            else
                if mode == "Jv"
                    HJv = zeros( Float64, m )
                else
                    HJv = zeros( Float64, n )
                end
                ic  = 0
                lder = length( pbm.conderlvl )
            end
        else
            if mode == "Jv"
                HJv = zeros( Float64, m )
            elseif mode == "Jtv"
                HJv = zeros( Float64, n )
            end
            ic  = 0
        end
    end

   # Check if pbm.A is not empty
   
    if !isempty(pbm.A)
        sA1, sA2 = size(pbm.A)
        has_A = true
    else
        has_A = false
    end
   
    # Evaluate the quadratic term, if any.

    if mode == "Hv" && !isempty(pbm.H)
        HJv += pbm.H * v
    end

    for iig in 1:length(glist)
        ig = glist[iig]
        
        #  Find the level of available derivatives for the group.

        if ( mode == "Jv" || mode == "Jtv" )
            if !isempty( pbm.conderlvl )
               if lder == 1 
                   derlvl = pbm.conderlvl[ 1 ]
               else
                   derlvl = pbm.conderlvl[ findfirst( isequal(ig), pbm.congrps ) ]
               end
            else
               derlvl = 2
            end

            #  Avoid computation for group ig if its first derivative is missing.
      
            if derlvl < 1 
                ic      = ic + 1
                HJv[ic] = NaN
                continue
            end
        end

    # Find the group's scaling.
    
        gsc = 1.0 # Default value
        if isdefined(pbm, :gscale)
            if ig <= length(pbm.gscale) && abs(pbm.gscale[ig]) > 1.0e-15
                gsc = pbm.gscale[ig]
            end
        end

        # Evaluate the linear term, if any.

        if isdefined(pbm, :gconst) && ig <= length(pbm.gconst)
            fin = -pbm.gconst[ig]
        else
            fin = 0.0
        end
        gin = zeros( Float64, n )
        if has_A && ig <= sA1
            fin += pbm.A[ig, 1:sA2]' * x[1:sA2]
            gin[1:sA2] = pbm.A[ig, 1:sA2]
        end
 
        Hin = spzeros( n, n)
 
        if debug 
            println("ig = $ig  fin(linear) = $fin")
        end
        
        #  Initialize the Hessian.

        Hinv = zeros( Float64, n )


        #  Loop on the group's elements.
        #
        #  The explicit sequential scalar assignment of gradient and Hessian components is
        #  necessary in this section because some elements may have repeated elemental variables,
        #  and only one occurence of such variables would be assigned by a vector assignment.
        
        if isdefined(pbm, :grelt)  && ig <= length(pbm.grelt) && !isempty(pbm.grelt[ig])
            for iiel in 1:length(pbm.grelt[ig])
                iel    = pbm.grelt[ig][iiel]
                irange = pbm.elvar[iel]
                efname = pbm.elftype[iel]
                if isdefined(pbm, :grelw) && ig <= length(pbm.grelw) && !isempty(pbm.grelw[ig])
                    has_weights = true
                    wiel = pbm.grelw[ig][iiel]
                else
                    has_weights = false
                end
 
                if mode == "Hv" || mode == "HIv"
                #  The group is an objective group.
                    fiel, giel, Hiel = pbm.call( efname, x[irange], iel, 3, pbm )
                    fin += has_weights ? wiel * fiel : fiel
                    for ir in 1:length(irange)
                        ii = irange[ir]
                        gin[ii] += has_weights ? wiel * giel[ir] : giel[ir]
                        for jr in 1:length(irange)
                           jj = irange[jr]
                           Hinv[ii] += has_weights ? wiel * Hiel[ir,jr]*v[jj] : Hiel[ir,jr]*v[jj]
                        end
                    end

                #  The group is an constraints group.

                elseif derlvl >= 1
                    fiel, giel = pbm.call( efname, x[irange], iel, 2, pbm )
                    fin += has_weights ? wiel * fiel : fiel
                    for ir in 1:length(irange)
                        ii = irange[ir]
                        gin[ii] += has_weights ? wiel * giel[ir] : giel[ir]
                    end
                end
            end
        end
 
        #  Evaluate the group function.
        #  1) the non-TRIVIAL case
        
        if isdefined(pbm, :grftype) && ig <= length(pbm.grftype) && isassigned(pbm.grftype,ig) && pbm.grftype[ig] != ""
            egname = pbm.grftype[ig]  # the group's ftype
        else 
            egname = "TRIVIAL"
        end
   
        if mode == "Hv"
            if (  egname == "TRIVIAL" )  
                HJv = HJv + Hinv / gsc
            else
                _, grada, Hessa = pbm.call( egname, fin, ig, 3, pbm )
                sgin = sparse(gin)
                HJv += ((Hessa * sgin) * dot( sgin , v ) + grada * Hinv) / gsc
            end
        elseif mode == "HIv"
            ic += 1
            if (  egname == "TRIVIAL" )  
                HJv = HJv + y[ic] * Hinv / gsc
            else
                _, grada, Hessa = pbm.call( egname, fin, ig, 3, pbm )
                sgin = sparse(gin)
                HJv += y[ic] * ((Hessa * sgin) * dot( sgin , v ) + grada * Hinv) / gsc
            end
        elseif mode == "Jv"
            ic += 1
            if derlvl >= 1
                if (  egname == "TRIVIAL" )
                    sgin    = sparse(gin)
                    HJv[ic] = dot( sgin , v ) / gsc
                else
                    _, grada = pbm.call( egname, fin, ig, 2, pbm )
                    sgin    = sparse(gin)
                    HJv[ic] = grada * dot( sgin , v ) / gsc
                end
            else
                HJv[ic] = NaN
            end
        elseif mode == "Jtv"
            ic += 1
            if derlvl >= 1
                if (  egname == "TRIVIAL" )
                    HJv  +=  gin * v[ic] / gsc
                else
                    _, grada = pbm.call( egname, fin, ig, 2, pbm )
                    HJv  += grada *  gin * v[ic] / gsc
                end
            else
                HJv[ic] = NaN
            end
         
        end

    end
    return HJv
end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#   Depending on mode:
#   mode = "Hv"  : evaluate the product of the objective's Hessian times v (glist = obj groups)
#   mode = "HIv" : evaluate the product of the constraints' Hessian  times v time the multiplier y
#                 (glist = cons groups)
#   mode = "Jv"  : evaluate the product of the constraints' Jacobian times v (glist = cons groups)
#   The vector y is unused (and unreferenced) for modes "Hv" and "Jv".
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function evalLx( gobjlist::Vector{Int}, gconlist::Vector{Int},
                 x::Vector{Float64}, y::Vector{Float64}, pbm::PBM, nargout::Int )

    # Handling the case of function evaluation
    
    if nargout == 1
        if length( gobjlist ) > 0 || !isempty( pbm.H )
            Lxy = evalgrsum( true, gobjlist, x, pbm, 1 )
        else
            Lxy = 0.0
        end
        if length( gconlist ) > 0
            c    = evalgrsum( false, gconlist, x, pbm, 1 )
            Lxy += only( y' * c )
        end
        return Float64(Lxy)

    # Handling the case of gradient evaluation
    
    elseif nargout == 2

        if length( gobjlist ) > 0 || !isempty( pbm.H )
            Lxy, Lgxy = evalgrsum( true, gobjlist, x, pbm, 2 )
        else
            Lxy  = 0.0
            Lgxy = zeros( Float64, length( x ) )
        end
        if length( gconlist ) > 0
            c, J  = evalgrsum( false, gconlist, x, pbm, 2 )
            Lxy  += only( y' * c )
            Lgxy += J' * y
        end
        return Float64(Lxy), Lgxy
        
    # Handling the case of Hessian evaluation
    
    elseif nargout == 3
    
        if length( gobjlist ) > 0 || !isempty( pbm.H )
            Lxy, Lgxy, LgHxy = evalgrsum( true, gobjlist, x, pbm, 3 )
        else
            n     = length( x )
            Lxy   = 0.0
            Lgxy  = zeros( Float64, n )
            LgHxy = spzeros( n, n )  # Use spzeros to create a sparse matrix
        end
        if length( gconlist ) > 0
            c, J, cHi = evalgrsum( false, gconlist, x, pbm, 3 )
            Lxy  += only( y'* c )
            Lgxy += J' * y
            for ig in 1:length( gconlist )
                LgHxy += y[ig] * cHi[ig]
            end
        end
        return Float64(Lxy), Lgxy, LgHxy
    end
    
end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#   Evaluate the product of the Lagrangian's Hessian times a vector v,
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function evalLHxyv( gobjlist::Vector{Int}, gconlist::Vector{Int},
                    x::Vector{Float64}, y::Vector{Float64}, v::Vector{Float64}, pbm::PBM )

    if length(gobjlist) > 0 || isdefined( pbm, :H)
        LHxyv = evalHJv( "Hv", gobjlist, x, v, Float64[], pbm )
    else
        n     = length( x )
        LHxyv = zeros( Float64, n )
    end
    if length( gconlist ) > 0
        for ig in 1: length( gconlist )
            LHxyv += evalHJv( "HIv", [gconlist[ig]], x, v, y, pbm )
        end
    end
    return LHxyv

end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#   This tool consider all problems in list_of_julia_problems (whose files are in
#   the ./julia_problems directory) and selects those whose SIF classification matches
#   that given by the input string classif. Matching is in the sense of  regular expressions
#   (regexp).
#   If varargin is empty (i.e. only classif is used as input argument), the function prints
#   the list of matching problems on standard output. Otherwise, the list is output in the
#   file whose name is a string passed as varargin{1} (Further components of varargin are
#   ignored).
#
#   If the input string is 'help'  or 'h', a message is printed on the standard output
#   describing the SIF classification scheme and an brief explanation of how to use the tool.
#
#   Thanks to Greta Malaspina (Firenze) for an inital implementation in Matlab.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function s2mpjlib_select( classif::String, args::Union{String,Any} = [] )

    if classif == "help" || classif == "h"

        println( "  " )
        println( " === The classification scheme ===" )
        println( "  " )
        println( " A problem is classified by a string of the form" )
        println( "    X-XXXXr-XX-n-m" )
        println( " The first character in the string identifies the problem collection" )
        println( " from which the problem is extracted. Possible values are" )
        println( "    C the CUTEst collection;" )
        println( "    S the SPARCO collection; and" )
        println( "    N none of the above." )
        println( " The character immediately following the first hyphen defines the type" )
        println( " of variables occurring in the problem. Its possible values are" )
        println( "    C the problem has continuous variables only;" )
        println( "    I the problem has integer variables only;" )
        println( "    B the problem has binary variables only; and" )
        println( "    M the problem has variables of different types." )
        println( " The second character after the first hyphen defines the type" )
        println( " of the problem''s objective function. Its possible values are" )
        println( "    N no objective function is defined;" )
        println( "    C the objective function is constant;" )
        println( "    L the objective function is linear;" )
        println( "    Q the objective function is quadratic;" )
        println( "    S the objective function is a sum of squares; and" )
        println( "    O the objective function is none of the above." )
        println( " The third character after the first hyphen defines the type of" )
        println( " constraints of the problem. Its possible values are" )
        println( "    U the problem is unconstrained;" )
        println( "    X the problem’s only constraints are fixed variables;" )
        println( "    B the problem’s only constraints are bounds on the variables;" )
        println( "    N the problem’s constraints represent the adjacency matrix of a (linear)" )
        println( "      network;" )
        println( "    L the problem’s constraints are linear;" )
        println( "    Q the problem’s constraints are quadratic; and" )
        println( "    O the problem’s constraints are more general than any of the above alone." )
        println( " The fourth character after the first hyphen indicates the smoothness of" )
        println( " the problem. There are two possible choices" )
        println( "    R the problem is regular, that is, its first and second derivatives " )
        println( "      exist and are continuous everywhere; or" )
        println( "    I the problem is irregular." )
        println( " The integer (r) which corresponds to the fourth character of the string is" )
        println( " the degree of the highest derivatives provided analytically within the problem" )
        println( " description. It is restricted to being one of the single characters O, 1, or 2." )
        println( " The character immediately following the second hyphen indicates the primary" )
        println( " origin and/or interest of the problem. Its possible values are" )
        println( "    A the problem is academic, that is, has been constructed specifically by" )
        println( "      researchers to test one or more algorithms;" )
        println( "    M the problem is part of a modeling exercise where the actual value of the" )
        println( "      solution is not used in a genuine practical application; and" )
        println( "    R the problem’s solution is (or has been) actually used in a real")
        println( "      application for purposes other than testing algorithms." )
        println( " The next character in the string indicates whether or not the problem" )
        println( " description contains explicit internal variables. There are two possible" )
        println( " values, namely," )
        println( "    Y the problem description contains explicit internal variables; or" )
        println( "    N the problem description does not contain any explicit internal variables." )
        println( " The symbol(s) between the third and fourth hyphen indicate the number of" )
        println( " variables in the problem. Possible values are" )
        println( "    V the number of variables in the problem can be chosen by the user; or" )
        println( "    n a positive integer giving the actual (fixed) number of problem variables." )
        println( " The symbol(s) after the fourth hyphen indicate the number of constraints" )
        println( " (other than fixed variables and bounds) in the problem. Note that fixed" )
        println( " variables are not considered as general constraints here. The two possible" )
        println( " values are" )
        println( "    V the number of constraints in the problem can be chosen by the user; or" )
        println( "    m a nonnegative integer giving the actual (fixed) number of constraints." )
        println( "  " )
        println( " === Using the problem selection tool ===" )
        println( "  " )
        println( " To use the selection tool, you should first open Julia in the parent of" )
        println( " the directory containing the Julia problem files, include the library by" )
        println( " issuing the command" )
        println( "    include( \"s2mpjlib.jl\" ) " )
        println( " and then call the selection function 's2mpjlib_select' itself." )
        println( " This function's first argument is a string specifying the class of problems" )
        println( " of interest.  This string is constructed by replacing by a dot each " )
        println( " character in the classification string for which all possible values are" )
        println( " acceptable (the dot is a wildcard character). For instance" )
        println( "    s2mpjlib_select( \"C-CSU..-..-2-0\" ) ")
        println( " lists all CUTEst unconstrained \"sum-of-squares\" problems in two continuous" )
        println( " variables, while " )
        println( "    s2mpjlib_select( \"C-C....-..-V-V\" ) " )
        println( " lists all CUTEst problems with variable number of continuous variables and" )
        println( " variable number of constraints." )
        println( " The classification strings \"unconstrained\", \"bound-constrained\", " )
        println( " \"fixed-variables\", \"general-constraints\", \"variable-n\" and " )
        println( " \"variable-m\" are also allowed." )
        println( " NOTE: any regular expression may be used as the first argument of " )
        println( "       s2mpj_select.problems to specify the problem class, so that, for" )
        println( "       instance, the previous selection can also be achieved by the command" )
        println( "            s2mpjlib_select( \"C-C.*V-V\" ) ")
        println( " Writing the list of selected problems to a file is obtained by specifying" )
        println( " the name of the file as a second argument of select, as in ")
        println( "    s2mpjlib_select( \"C-C....-..-V-V\", filename )" )

    else

        if classif == "unconstrained"
            classif = ".-..U.*"
        elseif classif == "bound-constrained"
            classif = ".-..B.*"
        elseif classif == "fixed-variables"
            classif = ".-..X.*"
        elseif classif == "general-constraints"
            classif = ".-..[LNQO].*"
        elseif classif == "variable-n"
            classif = ".-..B..-..-V-[V0-9]*"
        elseif classif == "variable-m"
            classif = ".-..B..-..-[V0-9]*-V"
        else
            lencl = length( classif )
            t = classif[12] == '.'
            if lencl > 11 && classif[12] == '.'
                oclassif = classif
                classif  = classif[1:11] * "[V0-9]*"
                if lencl> 12
                    classif = classif * oclassif[13:lencl]
                end
            end
            lenclm1 = length( classif ) - 1
            if classif[ lenclm1 ] == '.'
                classif = classif[1:lenclm1] * "[V0-9]*"
            end
        end
        filter = Regex( "classification = \"" * classif )

        list_of_problems = "./list_of_julia_problems"
        julia_problems   = "./julia_problems/"

        fid = ( length( args ) > 0 ) ? open( args, "w" ) : stdout
        
        allprobs = readlines( list_of_problems )
        
        for theprob in allprobs
            problem = string( julia_problems, theprob )
            if isfile( problem ) && occursin( filter, read( problem, String ) )
                println( fid, theprob ) 
            end
        end

        if fid != stdout
            close(fid)
        end
    end
    
end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

