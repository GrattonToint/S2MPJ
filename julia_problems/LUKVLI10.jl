function LUKVLI10(action::String,args::Union{Any}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : LUKVLI10
#    *********
# 
#    Source: Problem 5.10, the generalized Brown function with
#    Broyden tridiagonal constraints, due to L. Luksan and J. Vlcek,
#    "Sparse and partially separable test problems for 
#    unconstrained and equality constrained optimization",
#    Technical Report 767, Inst. Computer Science, Academy of Sciences
#    of the Czech Republic, 182 07 Prague, Czech Republic, 1999
# 
#    Equality constraints changed to inequalities
# 
#    SIF input: Nick Gould, April 2001
# 
#    classification = "C-COQR2-AY-V-V"
# 
#    some useful parameters, including N, the number of variables.
# 
#       Alternative values for the SIF file parameters:
# IE N                   100            $-PARAMETER
# IE N                   1000           $-PARAMETER
# IE N                   10000          $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 21 VI 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "LUKVLI10"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling LUKVLI10.")
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
            v_["N"] = Int64(10);  #  SIF file default value
        else
            v_["N"] = Int64(args[1]);
        end
# IE N                   100000         $-PARAMETER
        v_["1"] = 1
        v_["2"] = 2
        v_["N/2"] = trunc(Int,(v_["N"]/v_["2"]))
        v_["N-2"] = -2+v_["N"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        irA   = Int64[]
        icA   = Int64[]
        valA  = Float64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype = String[]
        for I = Int64(v_["1"]):Int64(v_["N/2"])
            ig,ig_,_ = s2mpj_ii("OBJ1"*string(I),ig_)
            arrset(gtype,ig,"<>")
            ig,ig_,_ = s2mpj_ii("OBJ2"*string(I),ig_)
            arrset(gtype,ig,"<>")
        end
        for K = Int64(v_["1"]):Int64(v_["N-2"])
            v_["K+1"] = 1+K
            v_["K+2"] = 2+K
            ig,ig_,_ = s2mpj_ii("C"*string(K),ig_)
            arrset(gtype,ig,"<=")
            arrset(pb.cnames,ig,"C"*string(K))
            push!(irA,ig)
            push!(icA,ix_["X"*string(Int64(v_["K+1"]))])
            push!(valA,Float64(3.0))
            push!(irA,ig)
            push!(icA,ix_["X"*string(K)])
            push!(valA,Float64(-1.0))
            push!(irA,ig)
            push!(icA,ix_["X"*string(Int64(v_["K+2"]))])
            push!(valA,Float64(-2.0))
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
        for K = Int64(v_["1"]):Int64(v_["N-2"])
            pbm.gconst[ig_["C"*string(K)]] = Float64(-1.0)
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        for I = Int64(v_["1"]):Int64(v_["2"]):Int64(v_["N"])
            pb.x0[ix_["X"*string(I)]] = Float64(-1.0)
        end
        for I = Int64(v_["2"]):Int64(v_["2"]):Int64(v_["N"])
            pb.x0[ix_["X"*string(I)]] = Float64(1.0)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSQR", iet_)
        loaset(elftv,it,1,"V")
        it,iet_,_ = s2mpj_ii( "eNASTY", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["N/2"])
            v_["2I"] = 2*I
            v_["2I-1"] = -1+v_["2I"]
            ename = "OBJ1"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eNASTY")
            arrset(ielftype,ie,iet_["eNASTY"])
            vname = "X"*string(Int64(v_["2I-1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["2I"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "OBJ2"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eNASTY")
            arrset(ielftype,ie,iet_["eNASTY"])
            vname = "X"*string(Int64(v_["2I"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["2I-1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        for K = Int64(v_["1"]):Int64(v_["N-2"])
            v_["K+1"] = 1+K
            ename = "C"*string(K)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
            vname = "X"*string(Int64(v_["K+1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N/2"])
            ig = ig_["OBJ1"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["OBJ1"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            ig = ig_["OBJ2"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["OBJ2"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        for K = Int64(v_["1"]):Int64(v_["N-2"])
            ig = ig_["C"*string(K)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["C"*string(K)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(-2.0))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN               3.53122E+02
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-COQR2-AY-V-V"
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
            g_[1] = 2.0*EV_[1]
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

    elseif action == "eNASTY"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        XX = EV_[1]*EV_[1]
        YYP1 = EV_[2]*EV_[2]+1.0
        F = XX^YYP1
        TLOGXX = 2.0*log(XX)
        FX = 2.0*YYP1/EV_[1]
        FY = TLOGXX*EV_[2]
        FXDOT = F*FX
        FYDOT = F*FY
        f_   = F
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = FXDOT
            g_[2] = FYDOT
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = FXDOT*FX-2.0*F*YYP1/XX
                H_[1,2] = FYDOT*FX+4.0*F*EV_[2]/EV_[1]
                H_[2,1] = H_[1,2]
                H_[2,2] = FYDOT*FY+F*TLOGXX
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

