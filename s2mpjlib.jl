#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                                S2MPJ library for Julia
#
#   Performs the runtime actions specific to S2MPJ, irrespective of the problem at hand.
#
#   Programming: S. Gratton and Ph. L. Toint (this version 11 VI 2024)
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
                            String[]   #cnames
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
    call::Function
    
    PBM(name::String) = new( name,                             #name
                             Int64[],                          #objgrps
                             Int64[],                          #congrps
                             spzeros(Float64, 100000, 100000), #A
                             Float64[],                        #gconst
                             spzeros(Float64, 100000, 100000), #H
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
                             [0, 0]                            #has_globs
                             )
end


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#   Computes the effective index of name in List, or add name to List if not in there already.
#   Return the index of name in the dictionary and new = 1 if the dictionary has been enlarged
#   or 0 otherwise.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s2mpj_ii( name, dict::Dict{String,Int} )
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

function s2mpj_nlx( name, dict::Dict{String,Int}, pb=nothing, getxnames=nothing,
                  xlowdef=nothing, xuppdef=nothing, x0def=nothing )

    iv, dict, newvar = s2mpj_ii( name, dict )

    if newvar == 1
        pb.n += 1
        if !isnothing( getxnames )
            arrset( pb.xnames, iv, name )
        end
        if hasproperty( pb, :xlower )
            arrset( pb.xlower, iv, isnothing( xlowdef ) ? 0.0 : xlowdef )
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
               push!(thelist, "" )  # Add elements to reach iv2
           end
        elseif isa( value, Int64)
           append!(thelist, zeros(Int,iv-length(thelist) ) )
        elseif isa( value, Float64)
           append!(thelist, zeros(Float64,iv-length(thelist) ) )
        else
           resize!(thelist, iv )
        end
    end
    thelist[iv] = value
    
end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#   Insert a value in a list of lists, possibly enlarging the list(s) to do so.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function loaset(list::Vector{Vector{T}}, iv1::Int, iv2::Int, value::T) where T
    while length(list) < iv1
        push!(list, Vector{T}())   # Add empty vectors to reach iv1
    end
    arrset( list[iv1], iv2, value )
#    if isa( value, String )
#        while length(list[iv1]) < iv2
#            push!(list[iv1], "" )  # Add elements to reach iv2
#        end
#    else
#        while length(list[iv1]) < iv2
#            push!(list[iv1], -1 )  # Add elements to reach iv2
#        end
#    end
#    list[iv1][iv2] = value
end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#   Perform the main high-level action requested from a S2MPJ problem file, using information computed
#   during the 'setup' call and passed to this function in the pbm struct.  The nonlinear functions
#   are evaluated at x = args[2]. See the documentation of S2MPJ for details.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s2mpj_eval( action, args... )

    pbm = args[1];

    #  Get the problem data and possible global parameters' values.
    
    if pbm.has_globs[1] >0 
        pbm = pbm.call( "e_globs", pbm )
    end
    if pbm.has_globs[2] >0 
        pbm = pbm.call( "g_globs", pbm )
    end
    
    #  Evaluations of the objective function and constraints
    
    if action == "fx"   
        return evalgrsum( true,  pbm.objgrps, args[2], pbm, 1 )
    elseif action == "fgx" 
        return evalgrsum( true,  pbm.objgrps, args[2], pbm, 2 )
    elseif action == "fgHx"
        return evalgrsum( true,  pbm.objgrps, args[2], pbm, 3 )
    elseif action == "cx"
        return evalgrsum( false, pbm.congrps, args[2], pbm, 1 )
    elseif action == "cJx"
        return evalgrsum( false, pbm.congrps, args[2], pbm, 2 )
    elseif action == "cJHx"
        return evalgrsum( false, pbm.congrps, args[2], pbm, 3 )
    elseif action == "cIx"
        return evalgrsum( false, pbm.congrps[args[3]], args[2], pbm, 1 )
    elseif action == "cIJx"
        return evalgrsum( false, pbm.congrps[args[3]], args[2], pbm, 2 )
    elseif action == "cIJHx"
        return evalgrsum( false, pbm.congrps[args[3]], args[2], pbm, 3 )
    elseif action == "fHxv"
        return evalHJv( "Hv", pbm.objgrps, args[2], args[3], [], pbm ); 
    elseif action == "cJxv"
        return evalHJv( "Jv", pbm.congrps, args[2], args[3], [], pbm );
    elseif action == "cIJxv"
        return evalHJv( "Jv", pbm.congrps[args[4]], args[2], args[3], [], pbm )
   
    # For the Lagrangian
    
    elseif action == "Lxy"
        return evalLx( pbm.objgrps, pbm.congrps, args[2], args[3], pbm, 1 );
    elseif action == "Lgxy"
        return evalLx( pbm.objgrps, pbm.congrps, args[2], args[3], pbm, 2 );
    elseif action == "LgHxy"
        return evalLx( pbm.objgrps, pbm.congrps, args[2], args[3], pbm, 3 );
    elseif action == "LIxy"
        return evalLx( pbm.objgrps, pbm.congrps[args[4]], args[2], args[3], pbm, 1 )
    elseif action == "LIgxy"
        return evalLx( pbm.objgrps, pbm.congrps[args[4]], args[2], args[3], pbm, 2 )
    elseif action == "LIgHxy"
        return evalLx( pbm.objgrps, pbm.congrps[args[4]], args[2], args[3], pbm, 3 )
    elseif action == "LHxyv"
        return evalLHxyv( pbm.objgrps, pbm.congrps, args[2], args[3], args[4], pbm )
    elseif action == "LIHxyv"
        return evalLHxyv( pbm.objgrps, pbm.congrps[args[5]], args[2], args[3], args[4], pbm )
    end
    
end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#   Evaluate the value of a sum of groups (and, if requested, that of of its gradient
#   and Hessian) at x, given the problem data available in the pbm struct.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function evalgrsum( isobj, glist, x, pbm, nargout )

   debug = false#D

#  Initializations

   n  = length( x );
   if ( isobj )
      fx = 0;
   else
      m  = length( glist );
      cx = zeros( m, 1 );
      ic = 0;
   end
   if ( nargout > 1 )
      if ( isobj )
         gx = zeros( n, 1 );
      else
         Jx = spzeros( m, n );
      end
      if ( nargout > 2 )
         if ( isobj )
            Hx = spzeros( n, n );
         else
            Hx = [];
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
           gin = zeros( n, 1 )
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
               if nargout == 1
                   fiel = pbm.call( efname, x[irange], iel, 1, pbm )
                   fin += has_weights ? wiel * fiel : fiel
               elseif nargout == 2
                   fiel, giel = pbm.call( efname, x[irange], iel, 2, pbm )
                   fin       += has_weights ? wiel * fiel : fiel
                   for ir in 1:length(irange)
                       ii       = irange[ir]
                       gin[ii] += has_weights ? wiel * giel[ir] : giel[ir]
                   end
               elseif nargout == 3
                   fiel, giel, Hiel = pbm.call( efname, x[irange], iel, 3, pbm )
#                   println( "irange = ", irange, "  x[irange] = ", x[irange], "  efname = ", efname, "  iel = ", iel, "  fiel = ", fiel, " wiel = ", wiel )#D
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
           
               # use evel(Meta.parse(.)) to replace the corresponding Matlab expression
               
               if nargout == 1
                   fa = pbm.call( egname, fin, ig, 1, pbm )
                   fx += fa / gsc
               elseif nargout == 2
                   fa, grada = pbm.call( egname, fin, ig, 2, pbm )
                   fx += fa / gsc
                   gx += grada[1] * gin / gsc
               elseif nargout == 3
#                   println( "fin = ", fin, " egname = ", egname, " gsc = ", gsc )#D
                   fa, grada, Hessa = pbm.call( egname, fin, ig, 3, pbm )
                   fx += fa / gsc
                   gx += grada[1] * gin / gsc
                   sgin = sparse(gin)
                   Hx += (Hessa[1] * sgin * sgin' + grada[1] * Hin) / gsc
               end
           else
               ic += 1
               if nargout == 1
                   fa = pbm.call( egname, fin, ig, 1, pbm )
                   cx[ic] = fa / gsc
               elseif nargout == 2
                   fa, grada = pbm.call( egname, fin, ig, 2, pbm )
                   cx[ic] = fa / gsc
                   Jx[ic, 1:n] = (grada[1] * gin') / gsc
               elseif nargout == 3
                   fa, grada, Hessa = pbm.call( egname, fin, ig, 3, pbm )
                   cx[ic] = fa / gsc
                   sgin = sparse(gin)
                   Jx[ic, 1:n] = (grada[1] * sgin') / gsc
                   push!(Hx, (Hessa[1] * sgin * sgin' + grada[1] * Hin) / gsc)
               end
           end
       else
       
           # TRIVIAL case
           
           if isobj
#               println( " objective fin = ", fin, " egname = TRIVIAL", " gsc = ", gsc )#D
               fx += fin / gsc
               if nargout >= 2
                   gx += gin / gsc
               end
               if nargout == 3
                   Hx += Hin / gsc
               end
           else
#               println( " constraint fin = ", fin, " egname = TRIVIAL", " gsc = ", gsc )#D
               ic += 1
               cx[ic] = fin / gsc
               if nargout >= 2
                   Jx[ic, 1:n] = gin' / gsc
               end
               if nargout == 3
                   push!( Hx, Hin / gsc)
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

function evalHJv( mode, glist, x, v, y, pbm )

    debug = false#D
   
    #  Initializations

    n  = length( x );
    if mode == "Hv"
        HJv = zeros( n, 1 );
    elseif mode == "HIv"
        HJv = zeros( n, 1 );
        ic  = 0;
    else
        HJv = zeros( length(glist), 1 );
        ic  = 0;
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
        HJv += pbm.H * v;
    end

    for iig in 1:length(glist)
        ig = glist[iig]
        
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
        gin = zeros(n, 1)
        if has_A && ig <= sA1
            fin += pbm.A[ig, 1:sA2]' * x[1:sA2]
            gin[1:sA2] = pbm.A[ig, 1:sA2]
        end
 
        Hin = spzeros(n, n)
 
        if debug 
            println("ig = $ig  fin(linear) = $fin")
        end
        #  Initialize the Hessian.

        Hinv = zeros( n, 1 );


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
                else
                #  The group is an constraints group.
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
                HJv = HJv + Hinv / gsc;
            else
                _, grada, Hessa = pbm.call( egname, fin, ig, 3, pbm )
                sgin = sparse(gin)
                HJv += ((Hessa * sgin) * dot( sgin , v ) + grada * Hinv) / gsc;
            end
        elseif mode == "HIv"
            ic += 1
            if (  egname == "TRIVIAL" )  
                HJv = HJv + y[ic] * Hinv / gsc;
            else
                _, grada, Hessa = pbm.call( egname, fin, ig, 3, pbm )
                sgin = sparse(gin)
                HJv += y[ic] * ((Hessa * sgin) * dot( sgin , v ) + grada * Hinv) / gsc;
            end
        else
            ic += 1
            if (  egname == "TRIVIAL" )
                sgin = sparse(gin)
                HJv[ic] =  dot( sgin , v ) / gsc;
            else
                _, grada = pbm.call( egname, fin, ig, 2, pbm )
                sgin = sparse(gin)
                HJv[ic]=  grada * dot( sgin , v ) / gsc;
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

function evalLx( gobjlist, gconlist, x, y, pbm, nargout )

    # Handling the case of function evaluation
    
    if nargout == 1
        if length( gobjlist ) > 0 || isdefined( pbm, :H )
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

        if length( gobjlist ) > 0 || isdefined( pbm, :H)
            Lxy, Lgxy = evalgrsum( true, gobjlist, x, pbm, 2 )
        else
            Lxy  = 0.0
            Lgxy = zeros( length( x ), 1 )
        end
        if length( gconlist ) > 0
            c, J  = evalgrsum( false, gconlist, x, pbm, 2 )
            Lxy  += only( y' * c )
            Lgxy += J' * y
        end
        return Float64(Lxy), Lgxy
        
    # Handling the case of Hessian evaluation
    
    elseif nargout == 3
    
        if length( gobjlist ) > 0 || isdefined( pbm, :H)
            Lxy, Lgxy, LgHxy = evalgrsum( true, gobjlist, x, pbm, 3 )
        else
            n     = length( x )
            Lxy   = 0.0
            Lgxy  = zeros( n, 1 )
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

function evalLHxyv( gobjlist, gconlist, x, y, v, pbm )

    if length(gobjlist) > 0 || isdefined( pbm, :H)
        LHxyv = evalHJv( "Hv", gobjlist, x, v, [], pbm )
    else
        n     = length( x )
        LHxyv = zeros( n, 1 )  
    end
    if length( gconlist ) > 0
        for ig in 1: length( gconlist )
            LHxyv += evalHJv( "HIv", gconlist[ig], x, v, y, pbm )
        end
    end
    return LHxyv

end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

