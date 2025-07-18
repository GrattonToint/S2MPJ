function ARGAUSS(action::String,args::Union{Any}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : ARGAUSS
#    *********
#    More''s gaussian problem in 3 variables
# 
#    Source: Problem 9 in
#    J.J. More', B.S. Garbow and K.E. Hillstrom,
#    "Testing Unconstrained Optimization Software",
#    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
# 
#    See also Buckley#28
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "C-CNOR2-AN-3-15"
# 
#    This function  is a nonlinear least squares with 15 groups.  Each
#    group has a nonlinear element.
# 
#    Number of groups
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 21 VI 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "ARGAUSS"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling ARGAUSS.")
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
        v_["15"] = 15
        v_["1"] = 1
        v_["8"] = 8.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        irA   = Int64[]
        icA   = Int64[]
        valA  = Float64[]
        iv,ix_,_ = s2mpj_ii("X1",ix_)
        arrset(pb.xnames,iv,"X1")
        iv,ix_,_ = s2mpj_ii("X2",ix_)
        arrset(pb.xnames,iv,"X2")
        iv,ix_,_ = s2mpj_ii("X3",ix_)
        arrset(pb.xnames,iv,"X3")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype = String[]
        for I = Int64(v_["1"]):Int64(v_["15"])
            ig,ig_,_ = s2mpj_ii("G"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"G"*string(I))
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
        pbm.gconst[ig_["G1"]] = Float64(0.0009)
        pbm.gconst[ig_["G2"]] = Float64(0.0044)
        pbm.gconst[ig_["G3"]] = Float64(0.0175)
        pbm.gconst[ig_["G4"]] = Float64(0.0540)
        pbm.gconst[ig_["G5"]] = Float64(0.1295)
        pbm.gconst[ig_["G6"]] = Float64(0.2420)
        pbm.gconst[ig_["G7"]] = Float64(0.3521)
        pbm.gconst[ig_["G8"]] = Float64(0.3989)
        pbm.gconst[ig_["G9"]] = Float64(0.3521)
        pbm.gconst[ig_["G10"]] = Float64(0.2420)
        pbm.gconst[ig_["G11"]] = Float64(0.1295)
        pbm.gconst[ig_["G12"]] = Float64(0.0540)
        pbm.gconst[ig_["G13"]] = Float64(0.0175)
        pbm.gconst[ig_["G14"]] = Float64(0.0044)
        pbm.gconst[ig_["G15"]] = Float64(0.0009)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        if haskey(ix_,"X1")
            pb.x0[ix_["X1"]] = Float64(0.4)
        else
            pb.y0[findfirst(x->x==ig_["X1"],pbm.congrps)] = Float64(0.4)
        end
        if haskey(ix_,"X2")
            pb.x0[ix_["X2"]] = Float64(1.0)
        else
            pb.y0[findfirst(x->x==ig_["X2"],pbm.congrps)] = Float64(1.0)
        end
        if haskey(ix_,"X3")
            pb.x0[ix_["X3"]] = Float64(0.0)
        else
            pb.y0[findfirst(x->x==ig_["X3"],pbm.congrps)] = Float64(0.0)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eGAUSS", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        loaset(elftv,it,3,"V3")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"T")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["15"])
            v_["REALI"] = Float64(I)
            v_["MREALI"] = -1.0*v_["REALI"]
            v_["8-I"] = v_["8"]+v_["MREALI"]
            v_["TI"] = 0.5*v_["8-I"]
            ename = "E"*string(I)
            ie,ie_,newelt = s2mpj_ii(ename,ie_)
            if newelt > 0
                arrset(pbm.elftype,ie,"eGAUSS")
                arrset(ielftype,ie,iet_["eGAUSS"])
            end
            vname = "X1"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="T",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["TI"]))
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["15"])
            ig = ig_["G"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-CNOR2-AN-3-15"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eGAUSS"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        TMV3 = pbm.elpar[iel_][1]-EV_[3]
        TMV3SQ = -0.5e0*TMV3*TMV3
        EXPA = exp(EV_[2]*TMV3SQ)
        FVAL = EV_[1]*EXPA
        V2TMV3 = EV_[2]*TMV3
        f_   = FVAL
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EXPA
            g_[2] = TMV3SQ*FVAL
            g_[3] = V2TMV3*FVAL
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,2] = TMV3SQ*EXPA
                H_[2,1] = H_[1,2]
                H_[1,3] = V2TMV3*EXPA
                H_[3,1] = H_[1,3]
                H_[2,2] = TMV3SQ*TMV3SQ*FVAL
                H_[2,3] = (TMV3SQ*V2TMV3+TMV3)*FVAL
                H_[3,2] = H_[2,3]
                H_[3,3] = (V2TMV3*V2TMV3-EV_[2])*FVAL
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

