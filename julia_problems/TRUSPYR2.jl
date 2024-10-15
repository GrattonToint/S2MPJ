function TRUSPYR2(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    This is a structural optimization problem.
#    The problem is to minimize the weight of a given
#    8-bar truss structure formed as a pyramid for a given external load.
#    There are upper bounds on the normal stresses in the
#    bars and lower bounds on the cross-sectional areas of the bars.
# 
#    Source:
#    K. Svanberg, 
#    "Local and global optima", 
#    Proceedings of the NATO/DFG ASI on Optimization of large structural
#    systems, 
#    G. I. N. Rozvany, ed., Kluwer, 1993, pp. 579-588.
# 
#    SIF input: A. Forsgren, Royal Institute of Technology, December 1993.
# 
#    classification = "C-LQR2-MN-11-11"
# 
#    Number of bars
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "TRUSPYR2"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["NBAR"] = 8
        v_["NDIM"] = 3
        v_["1"] = 1
        v_["2"] = 2
        v_["NBAR/2"] = trunc(Int,(v_["NBAR"]/v_["2"]))
        v_["8.0"] = 8.0
        v_["SQRT17"] = sqrt(17.0)
        v_["SQRT18"] = sqrt(18.0)
        v_["P1"] = 40.0
        v_["P2"] = 20.0
        v_["P3"] = 200.0
        for J = Int64(v_["1"]):Int64(v_["NBAR/2"])
            v_["L"*string(J)] = v_["SQRT17"]/v_["8.0"]
            v_["J+4"] = J+v_["NBAR/2"]
            v_["L"*string(Int64(v_["J+4"]))] = v_["SQRT18"]/v_["8.0"]
        end
        v_["E"] = 21.0
        v_["R1,1"] = 0.250
        v_["R2,1"] = 0.250
        v_["R3,1"] = 0.375
        v_["R1,2"] = 0.250
        v_["R2,2"] = -0.250
        v_["R3,2"] = 0.375
        v_["R1,3"] = -0.250
        v_["R2,3"] = -0.250
        v_["R3,3"] = 0.375
        v_["R1,4"] = -0.250
        v_["R2,4"] = 0.250
        v_["R3,4"] = 0.375
        v_["R1,5"] = 0.375
        v_["R2,5"] = 0.000
        v_["R3,5"] = 0.375
        v_["R1,6"] = 0.000
        v_["R2,6"] = -0.375
        v_["R3,6"] = 0.375
        v_["R1,7"] = -0.375
        v_["R2,7"] = 0.000
        v_["R3,7"] = 0.375
        v_["R1,8"] = 0.000
        v_["R2,8"] = 0.375
        v_["R3,8"] = 0.375
        for J = Int64(v_["1"]):Int64(v_["NBAR"])
            v_["L2"*string(J)] = v_["L"*string(J)]*v_["L"*string(J)]
            v_["L3"*string(J)] = v_["L2"*string(J)]*v_["L"*string(J)]
            v_["GAMMA"*string(J)] = v_["E"]/v_["L3"*string(J)]
            v_["DL2"*string(J)] = v_["L2"*string(J)]/v_["E"]
            v_["W"*string(J)] = 0.78*v_["L"*string(J)]
            v_["STRUP"*string(J)] = 10.0*v_["DL2"*string(J)]
            for I = Int64(v_["1"]):Int64(v_["NDIM"])
                v_["RG"*string(I)*","*string(J)]  = (
                      v_["GAMMA"*string(J)]*v_["R"*string(I)*","*string(J)])
                for K = Int64(v_["1"]):Int64(v_["NDIM"])
                    v_["RR"*string(I)*","*string(J)*","*string(K)]  = (
                          v_["RG"*string(I)*","*string(J)]*v_["R"*string(K)*","*string(J)])
                end
            end
        end
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for J = Int64(v_["1"]):Int64(v_["NBAR"])
            iv,ix_,_ = s2mpj_ii("XAREA"*string(J),ix_)
            arrset(pb.xnames,iv,"XAREA"*string(J))
        end
        for I = Int64(v_["1"]):Int64(v_["NDIM"])
            iv,ix_,_ = s2mpj_ii("DISPL"*string(I),ix_)
            arrset(pb.xnames,iv,"DISPL"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for J = Int64(v_["1"]):Int64(v_["NBAR"])
            ig,ig_,_ = s2mpj_ii("WEIGHT",ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["XAREA"*string(J)]
            pbm.A[ig,iv] += Float64(v_["W"*string(J)])
        end
        for K = Int64(v_["1"]):Int64(v_["NDIM"])
            ig,ig_,_ = s2mpj_ii("EQUIL"*string(K),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"EQUIL"*string(K))
        end
        for I = Int64(v_["1"]):Int64(v_["NDIM"])
            for J = Int64(v_["1"]):Int64(v_["NBAR"])
                ig,ig_,_ = s2mpj_ii("STRES"*string(J),ig_)
                arrset(gtype,ig,"<=")
                arrset(pb.cnames,ig,"STRES"*string(J))
                iv = ix_["DISPL"*string(I)]
                pbm.A[ig,iv] += Float64(v_["R"*string(I)*","*string(J)])
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
        for K = Int64(v_["1"]):Int64(v_["NDIM"])
            pbm.gconst[ig_["EQUIL"*string(K)]] = Float64(v_["P"*string(K)])
        end
        for J = Int64(v_["1"]):Int64(v_["NBAR"])
            pbm.gconst[ig_["STRES"*string(J)]] = Float64(v_["STRUP"*string(J)])
        end
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        for J = Int64(v_["1"]):Int64(v_["NBAR"])
            pb.xlower[ix_["XAREA"*string(J)]] = 1.0
        end
        for I = Int64(v_["1"]):Int64(v_["NDIM"])
            pb.xlower[ix_["DISPL"*string(I)]] = -Inf
            pb.xupper[ix_["DISPL"*string(I)]] = +Inf
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "en2PR", iet_)
        loaset(elftv,it,1,"U")
        loaset(elftv,it,2,"X")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["NDIM"])
            for J = Int64(v_["1"]):Int64(v_["NBAR"])
                ename = "UX"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"en2PR")
                arrset(ielftype,ie,iet_["en2PR"])
                vname = "DISPL"*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="U",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "XAREA"*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="X",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["NDIM"])
            for J = Int64(v_["1"]):Int64(v_["NBAR"])
                for K = Int64(v_["1"]):Int64(v_["NDIM"])
                    ig = ig_["EQUIL"*string(K)]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["UX"*string(I)*","*string(J)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RR"*string(I)*","*string(J)*","*string(K)]))
                end
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Objective function value corresponding to the local minimizer above
        pb.objlower = 1.2287408808
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-LQR2-MN-11-11"
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

    elseif action == "en2PR"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[2]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[2]
            g_[2] = EV_[1]
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = 1.0
                H_[2,1] = H_[1,2]
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

