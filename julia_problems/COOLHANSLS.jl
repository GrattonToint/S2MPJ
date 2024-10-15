function COOLHANSLS(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : COOLHANSLS
#    *********
# 
#    A problem arising from the analysis of a Cooley-Hansen economy with
#    loglinear approximation.  The problem is to solve the matrix equation
#                  A * X * X + B * X + C = 0
#    where A, B and C are known N times N matrices and X an unknown matrix
#    of matching dimension.  The instance considered here has N = 3.
# 
#    Source:
#    S. Ceria, private communication, 1995.
# 
#    SIF input: Ph. Toint, Feb 1995.
#    Least-squares version of COOLHANS.SIF, Nick Gould, Jan 2020.
# 
#    classification = "C-SUR2-RN-9-0"
# 
#    order of the matrix equation
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "COOLHANSLS"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["N"] = 3
        v_["A1,1"] = 0.0
        v_["A2,1"] = 0.13725e-6
        v_["A3,1"] = 0.0
        v_["A1,2"] = 0.0
        v_["A2,2"] = 937.62
        v_["A3,2"] = 0.0
        v_["A1,3"] = 0.0
        v_["A2,3"] = -42.207
        v_["A3,3"] = 0.0
        v_["B1,1"] = 0.0060893
        v_["B2,1"] = 0.13880e-6
        v_["B3,1"] = -0.13877e-6
        v_["B1,2"] = -44.292
        v_["B2,2"] = -1886.0
        v_["B3,2"] = 42.362
        v_["B1,3"] = 2.0011
        v_["B2,3"] = 42.362
        v_["B3,3"] = -2.0705
        v_["C1,1"] = 0.0
        v_["C2,1"] = 0.0
        v_["C3,1"] = 0.0
        v_["C1,2"] = 44.792
        v_["C2,2"] = 948.21
        v_["C3,2"] = -42.684
        v_["C1,3"] = 0.0
        v_["C2,3"] = 0.0
        v_["C3,3"] = 0.0
        v_["1"] = 1
        v_["2"] = 2
        v_["3"] = 3
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            for J = Int64(v_["1"]):Int64(v_["N"])
                iv,ix_,_ = s2mpj_ii("X"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"X"*string(I)*","*string(J))
            end
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for K = Int64(v_["1"]):Int64(v_["N"])
            for L = Int64(v_["1"]):Int64(v_["N"])
                for M = Int64(v_["1"]):Int64(v_["N"])
                    ig,ig_,_ = s2mpj_ii("G"*string(K)*","*string(L),ig_)
                    arrset(gtype,ig,"<>")
                    iv = ix_["X"*string(M)*","*string(L)]
                    pbm.A[ig,iv] += Float64(v_["B"*string(K)*","*string(M)])
                end
            end
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        for K = Int64(v_["1"]):Int64(v_["N"])
            for L = Int64(v_["1"]):Int64(v_["N"])
                v_["-C"] = -1.0*v_["C"*string(K)*","*string(L)]
                pbm.gconst[ig_["G"*string(K)*","*string(L)]] = Float64(v_["-C"])
            end
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "en2PR", iet_)
        loaset(elftv,it,1,"XX")
        loaset(elftv,it,2,"YY")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for K = Int64(v_["1"]):Int64(v_["N"])
            for L = Int64(v_["1"]):Int64(v_["N"])
                for M = Int64(v_["1"]):Int64(v_["N"])
                    ename = "E"*string(K)*","*string(M)*","*string(L)
                    ie,ie_,_  = s2mpj_ii(ename,ie_)
                    arrset(pbm.elftype,ie,"en2PR")
                    arrset(ielftype,ie,iet_["en2PR"])
                    vname = "X"*string(K)*","*string(M)
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    vname = "X"*string(M)*","*string(L)
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="YY",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                end
            end
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gL2",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for ig in 1:ngrp
            arrset(pbm.grftype,ig,"gL2")
        end
        for L = Int64(v_["1"]):Int64(v_["N"])
            for P = Int64(v_["1"]):Int64(v_["N"])
                for M = Int64(v_["1"]):Int64(v_["N"])
                    ig = ig_["G"*string(Int64(v_["1"]))*","*string(L)]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["E"*string(P)*","*string(M)*","*string(L)])
                    loaset(pbm.grelw,ig,posel,Float64(v_["A"*string(Int64(v_["1"]))*","*string(P)]))
                end
            end
        end
        for L = Int64(v_["1"]):Int64(v_["N"])
            for P = Int64(v_["1"]):Int64(v_["N"])
                for M = Int64(v_["1"]):Int64(v_["N"])
                    ig = ig_["G"*string(Int64(v_["2"]))*","*string(L)]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["E"*string(P)*","*string(M)*","*string(L)])
                    loaset(pbm.grelw,ig,posel,Float64(v_["A"*string(Int64(v_["2"]))*","*string(P)]))
                end
            end
        end
        for L = Int64(v_["1"]):Int64(v_["N"])
            for P = Int64(v_["1"]):Int64(v_["N"])
                for M = Int64(v_["1"]):Int64(v_["N"])
                    ig = ig_["G"*string(Int64(v_["3"]))*","*string(L)]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["E"*string(P)*","*string(M)*","*string(L)])
                    loaset(pbm.grelw,ig,posel,Float64(v_["A"*string(Int64(v_["3"]))*","*string(P)]))
                end
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN                0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-SUR2-RN-9-0"
        pb.x0          = zeros(Float64,pb.n)
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
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

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    elseif action == "gL2"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_= GVAR_*GVAR_
        if nargout>1
            g_ = GVAR_+GVAR_
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = 2.0
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

