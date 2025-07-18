function LIPPERT1(action::String,args::Union{Any}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : LIPPERT1
#    *********
# 
#    A discrete approximation to a continuum optimal flow problem
#    in the unit square. The continuum problem requires that the
#    divergence of a given flow should be given everywhere in the
#    region of interest, with the restriction that the capacity of
#    the flow is bounded. The aim is then to maximize the given flow.
# 
#    The discrete problem (primal formulation 1) in the unit square is to 
#      maximize   t
#      subject to dx( u_ij - ui-1j ) + dx( v_ij - vij-1 ) = t s_ij
#                 u_ij^2 + v_ij^2 <= 1
#                 u_i-1j^2 + v_ij^2 <= 1
#                 u_ij^2 + v_ij-1^2 <= 1
#                 u_i-1j^2 + v_ij-1^2 <= 1
#      where 1 <= i <= nx, 1 <= j <= ny
#      and        t >= 0
# 
#    Source: R. A. Lippert
#      "Discrete approximations to continuum optimal flow problems"
#      Tech. Report, Dept of Maths, M.I.T., 2006
#    following a suggestion by Gil Strang
# 
#    SIF input: Nick Gould, September 2006
# 
#    classification = "C-CLQR2-MN-V-V"
# 
#    Number of nodes in x direction
# 
#       Alternative values for the SIF file parameters:
# IE NX                  2              $-PARAMETER
# IE NX                  3              $-PARAMETER
# IE NX                  10             $-PARAMETER
# IE NX                  40             $-PARAMETER
# IE NX                  100            $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 21 VI 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "LIPPERT1"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling LIPPERT1.")
    end

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        if nargin<1
            v_["NX"] = Int64(3);  #  SIF file default value
        else
            v_["NX"] = Int64(args[1]);
        end
# IE NY                  2              $-PARAMETER
# IE NY                  3              $-PARAMETER
# IE NY                  10             $-PARAMETER 
# IE NY                  40             $-PARAMETER
# IE NY                  100            $-PARAMETER
        if nargin<2
            v_["NY"] = Int64(10);  #  SIF file default value
        else
            v_["NY"] = Int64(args[2]);
        end
        v_["X+"] = 1+v_["NX"]
        v_["X-"] = -1+v_["NX"]
        v_["Y+"] = 1+v_["NY"]
        v_["Y-"] = -1+v_["NY"]
        v_["1"] = 1
        v_["0"] = 0
        v_["ONE"] = 1.0
        v_["-ONE"] = -1.0
        v_["S"] = 1.0
        v_["-S"] = v_["S"]*v_["-ONE"]
        v_["RX"] = Float64(v_["NX"])
        v_["DX"] = v_["ONE"]/v_["RX"]
        v_["-DX"] = v_["-ONE"]/v_["RX"]
        v_["RY"] = Float64(v_["NY"])
        v_["DY"] = v_["ONE"]/v_["RY"]
        v_["-DY"] = v_["-ONE"]/v_["RY"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        irA   = Int64[]
        icA   = Int64[]
        valA  = Float64[]
        iv,ix_,_ = s2mpj_ii("T",ix_)
        arrset(pb.xnames,iv,"T")
        for I = Int64(v_["0"]):Int64(v_["NX"])
            for J = Int64(v_["1"]):Int64(v_["NY"])
                iv,ix_,_ = s2mpj_ii("U"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"U"*string(I)*","*string(J))
            end
        end
        for I = Int64(v_["1"]):Int64(v_["NX"])
            for J = Int64(v_["0"]):Int64(v_["NY"])
                iv,ix_,_ = s2mpj_ii("V"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"V"*string(I)*","*string(J))
            end
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["T"])
        push!(valA,Float64(-1.0))
        for I = Int64(v_["1"]):Int64(v_["NX"])
            v_["I-1"] = -1+I
            for J = Int64(v_["1"]):Int64(v_["NY"])
                v_["J-1"] = -1+J
                ig,ig_,_ = s2mpj_ii("O"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"O"*string(I)*","*string(J))
                arrset(pbm.gscale,ig,Float64(v_["DX"]))
                push!(irA,ig)
                push!(icA,ix_["U"*string(I)*","*string(J)])
                push!(valA,Float64(v_["DX"]))
                push!(irA,ig)
                push!(icA,ix_["U"*string(Int64(v_["I-1"]))*","*string(J)])
                push!(valA,Float64(v_["-DX"]))
                push!(irA,ig)
                push!(icA,ix_["V"*string(I)*","*string(J)])
                push!(valA,Float64(v_["DY"]))
                push!(irA,ig)
                push!(icA,ix_["V"*string(I)*","*string(Int64(v_["J-1"]))])
                push!(valA,Float64(v_["-DY"]))
                push!(irA,ig)
                push!(icA,ix_["T"])
                push!(valA,Float64(v_["-S"]))
            end
        end
        for I = Int64(v_["1"]):Int64(v_["NX"])
            for J = Int64(v_["1"]):Int64(v_["NY"])
                ig,ig_,_ = s2mpj_ii("A"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<=")
                arrset(pb.cnames,ig,"A"*string(I)*","*string(J))
                ig,ig_,_ = s2mpj_ii("B"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<=")
                arrset(pb.cnames,ig,"B"*string(I)*","*string(J))
                ig,ig_,_ = s2mpj_ii("C"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<=")
                arrset(pb.cnames,ig,"C"*string(I)*","*string(J))
                ig,ig_,_ = s2mpj_ii("D"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<=")
                arrset(pb.cnames,ig,"D"*string(I)*","*string(J))
            end
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        legrps = findall(x->x=="<=",gtype)
        eqgrps = findall(x->x=="==",gtype)
        gegrps = findall(x->x==">=",gtype)
        pb.nle = length(legrps)
        pb.neq = length(eqgrps)
        pb.nge = length(gegrps)
        pb.m   = pb.nle+pb.neq+pb.nge
        pbm.congrps = [[legrps;eqgrps];gegrps]
        pb.nob = ngrp-pb.m
        pbm.objgrps = findall(x->x=="<>",gtype)
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        for I = Int64(v_["1"]):Int64(v_["NX"])
            for J = Int64(v_["1"]):Int64(v_["NY"])
                pbm.gconst[ig_["A"*string(I)*","*string(J)]] = Float64(1.0)
                pbm.gconst[ig_["B"*string(I)*","*string(J)]] = Float64(1.0)
                pbm.gconst[ig_["C"*string(I)*","*string(J)]] = Float64(1.0)
                pbm.gconst[ig_["D"*string(I)*","*string(J)]] = Float64(1.0)
            end
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower[ix_["T"]] = 0.01
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSQR", iet_)
        loaset(elftv,it,1,"ALPHA")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["0"]):Int64(v_["NX"])
            for J = Int64(v_["1"]):Int64(v_["NY"])
                ename = "P"*string(I)*","*string(J)
                ie,ie_,newelt = s2mpj_ii(ename,ie_)
                if newelt > 0
                    arrset(pbm.elftype,ie,"eSQR")
                    arrset(ielftype,ie,iet_["eSQR"])
                end
                vname = "U"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="ALPHA",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        for I = Int64(v_["1"]):Int64(v_["NX"])
            for J = Int64(v_["0"]):Int64(v_["NY"])
                ename = "Q"*string(I)*","*string(J)
                ie,ie_,newelt = s2mpj_ii(ename,ie_)
                if newelt > 0
                    arrset(pbm.elftype,ie,"eSQR")
                    arrset(ielftype,ie,iet_["eSQR"])
                end
                vname = "V"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="ALPHA",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["NX"])
            v_["I-1"] = -1+I
            for J = Int64(v_["1"]):Int64(v_["NY"])
                v_["J-1"] = -1+J
                ig = ig_["A"*string(I)*","*string(J)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["P"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["Q"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
                ig = ig_["B"*string(I)*","*string(J)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["P"*string(Int64(v_["I-1"]))*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["Q"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
                ig = ig_["C"*string(I)*","*string(J)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["P"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["Q"*string(I)*","*string(Int64(v_["J-1"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
                ig = ig_["D"*string(I)*","*string(J)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["P"*string(Int64(v_["I-1"]))*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["Q"*string(I)*","*string(Int64(v_["J-1"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               -3.77245385
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-CLQR2-MN-V-V"
        pb.x0          = zeros(Float64,pb.n)
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eSQR"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[1]+EV_[1]
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 2.0
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    #%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    elseif action in  ["fx","fgx","fgHx","cx","cJx","cJHx","cIx","cIJx","cIJHx","cIJxv","fHxv",
                       "cJxv","cJtxv","cIJtxv","Lxy","Lgxy","LgHxy","LIxy","LIgxy","LIgHxy",
                       "LHxyv","LIHxyv"]

        pbm = args[1]
        if pbm.name == name
            pbm.has_globs = [0,0]
            return s2mpj_eval(action,args...)
        else
            println("ERROR: please run "*name*" with action = setup")
            return ntuple(i->undef,args[end])
        end

    else
        println("ERROR: action "*action*" unavailable for problem "*name*".jl")
        return ntuple(i->undef,args[end])
    end

end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

