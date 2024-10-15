function BROWNAL(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : BROWNAL
#    *********
#    Brown almost linear least squares problem.
#    This problem is a sum of n least-squares groups, the last one of
#    which has a nonlinear element.
#    It Hessian matrix is dense.
# 
#    Source: Problem 27 in
#    J.J. More', B.S. Garbow and K.E. Hillstrom,
#    "Testing Unconstrained Optimization Software",
#    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
# 
#    See also Buckley#79
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "C-SUR2-AN-V-0"
# 
#    N is the number of free variables (variable).
# 
#       Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER     original value
# IE N                   100            $-PARAMETER
# IE N                   200            $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "BROWNAL"

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
# IE N                   1000           $-PARAMETER
        v_["1"] = 1
        v_["N-1"] = -1+v_["N"]
        v_["N+1"] = 1+v_["N"]
        v_["RN+1"] = Float64(v_["N+1"])
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            v_["I-1"] = -1+I
            v_["I+1"] = 1+I
            for J = Int64(v_["1"]):Int64(v_["I-1"])
                ig,ig_,_ = s2mpj_ii("G"*string(I),ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["X"*string(J)]
                pbm.A[ig,iv] += Float64(1.0)
            end
            ig,ig_,_ = s2mpj_ii("G"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(2.0)
            for J = Int64(v_["I+1"]):Int64(v_["N"])
                ig,ig_,_ = s2mpj_ii("G"*string(I),ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["X"*string(J)]
                pbm.A[ig,iv] += Float64(1.0)
            end
        end
        ig,ig_,_ = s2mpj_ii("G"*string(Int64(v_["N"])),ig_)
        arrset(gtype,ig,"<>")
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            pbm.gconst[ig_["G"*string(I)]] = Float64(v_["RN+1"])
        end
        pbm.gconst[ig_["G"*string(Int64(v_["N"]))]] = Float64(1.0)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        for I = Int64(v_["1"]):Int64(v_["N"])
            pb.x0[ix_["X"*string(I)]] = Float64(0.5)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "ePROD", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        loaset(elftv,it,3,"V3")
        loaset(elftv,it,4,"V4")
        loaset(elftv,it,5,"V5")
        loaset(elftv,it,6,"V6")
        loaset(elftv,it,7,"V7")
        loaset(elftv,it,8,"V8")
        loaset(elftv,it,9,"V9")
        loaset(elftv,it,10,"V10")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "E"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD")
        arrset(ielftype,ie,iet_["ePROD"])
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
        vname = "X4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V5",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X6"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V6",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X7"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V7",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X8"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V8",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X9"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V9",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X10"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V10",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
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
        ig = ig_["G"*string(Int64(v_["N"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"])
        loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN               0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-SUR2-AN-V-0"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "ePROD"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        V12 = EV_[1]*EV_[2]
        V34 = EV_[3]*EV_[4]
        V56 = EV_[5]*EV_[6]
        V78 = EV_[7]*EV_[8]
        V910 = EV_[9]*EV_[10]
        f_   = V12*V34*V56*V78*V910
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[2]*V34*V56*V78*V910
            g_[2] = EV_[1]*V34*V56*V78*V910
            g_[3] = V12*EV_[4]*V56*V78*V910
            g_[4] = V12*EV_[3]*V56*V78*V910
            g_[5] = V12*V34*EV_[6]*V78*V910
            g_[6] = V12*V34*EV_[5]*V78*V910
            g_[7] = V12*V34*V56*EV_[8]*V910
            g_[8] = V12*V34*V56*EV_[7]*V910
            g_[9] = V12*V34*V56*V78*EV_[10]
            g_[10] = V12*V34*V56*V78*EV_[9]
            if nargout>2
                H_ = zeros(Float64,10,10)
                H_[1,2] = V34*V56*V78*V910
                H_[2,1] = H_[1,2]
                H_[1,3] = EV_[2]*EV_[4]*V56*V78*V910
                H_[3,1] = H_[1,3]
                H_[1,4] = EV_[2]*EV_[3]*V56*V78*V910
                H_[4,1] = H_[1,4]
                H_[1,5] = EV_[2]*V34*EV_[6]*V78*V910
                H_[5,1] = H_[1,5]
                H_[1,6] = EV_[2]*V34*EV_[5]*V78*V910
                H_[6,1] = H_[1,6]
                H_[1,7] = EV_[2]*V34*V56*EV_[8]*V910
                H_[7,1] = H_[1,7]
                H_[1,8] = EV_[2]*V34*V56*EV_[7]*V910
                H_[8,1] = H_[1,8]
                H_[1,9] = EV_[2]*V34*V56*V78*EV_[10]
                H_[9,1] = H_[1,9]
                H_[1,10] = EV_[2]*V34*V56*V78*EV_[9]
                H_[10,1] = H_[1,10]
                H_[2,3] = EV_[1]*EV_[4]*V56*V78*V910
                H_[3,2] = H_[2,3]
                H_[2,4] = EV_[1]*EV_[3]*V56*V78*V910
                H_[4,2] = H_[2,4]
                H_[2,5] = EV_[1]*V34*EV_[6]*V78*V910
                H_[5,2] = H_[2,5]
                H_[2,6] = EV_[1]*V34*EV_[5]*V78*V910
                H_[6,2] = H_[2,6]
                H_[2,7] = EV_[1]*V34*V56*EV_[8]*V910
                H_[7,2] = H_[2,7]
                H_[2,8] = EV_[1]*V34*V56*EV_[7]*V910
                H_[8,2] = H_[2,8]
                H_[2,9] = EV_[1]*V34*V56*V78*EV_[10]
                H_[9,2] = H_[2,9]
                H_[2,10] = EV_[1]*V34*V56*V78*EV_[9]
                H_[10,2] = H_[2,10]
                H_[3,4] = V12*V56*V78*V910
                H_[4,3] = H_[3,4]
                H_[3,5] = V12*EV_[4]*EV_[6]*V78*V910
                H_[5,3] = H_[3,5]
                H_[3,6] = V12*EV_[4]*EV_[5]*V78*V910
                H_[6,3] = H_[3,6]
                H_[3,7] = V12*EV_[4]*V56*EV_[8]*V910
                H_[7,3] = H_[3,7]
                H_[3,8] = V12*EV_[4]*V56*EV_[7]*V910
                H_[8,3] = H_[3,8]
                H_[3,9] = V12*EV_[4]*V56*V78*EV_[10]
                H_[9,3] = H_[3,9]
                H_[3,10] = V12*EV_[4]*V56*V78*EV_[9]
                H_[10,3] = H_[3,10]
                H_[4,5] = V12*EV_[3]*EV_[6]*V78*V910
                H_[5,4] = H_[4,5]
                H_[4,6] = V12*EV_[3]*EV_[5]*V78*V910
                H_[6,4] = H_[4,6]
                H_[4,7] = V12*EV_[3]*V56*EV_[8]*V910
                H_[7,4] = H_[4,7]
                H_[4,8] = V12*EV_[3]*V56*EV_[7]*V910
                H_[8,4] = H_[4,8]
                H_[4,9] = V12*EV_[3]*V56*V78*EV_[10]
                H_[9,4] = H_[4,9]
                H_[4,10] = V12*EV_[3]*V56*V78*EV_[9]
                H_[10,4] = H_[4,10]
                H_[5,6] = V12*V34*V78*V910
                H_[6,5] = H_[5,6]
                H_[5,7] = V12*V34*EV_[6]*EV_[8]*V910
                H_[7,5] = H_[5,7]
                H_[5,8] = V12*V34*EV_[6]*EV_[7]*V910
                H_[8,5] = H_[5,8]
                H_[5,9] = V12*V34*EV_[6]*V78*EV_[10]
                H_[9,5] = H_[5,9]
                H_[5,10] = V12*V34*EV_[6]*V78*EV_[9]
                H_[10,5] = H_[5,10]
                H_[6,7] = V12*V34*EV_[5]*EV_[8]*V910
                H_[7,6] = H_[6,7]
                H_[6,8] = V12*V34*EV_[5]*EV_[7]*V910
                H_[8,6] = H_[6,8]
                H_[6,9] = V12*V34*EV_[5]*V78*EV_[10]
                H_[9,6] = H_[6,9]
                H_[6,10] = V12*V34*EV_[5]*V78*EV_[9]
                H_[10,6] = H_[6,10]
                H_[7,8] = V12*V34*V56*V910
                H_[8,7] = H_[7,8]
                H_[7,9] = V12*V34*V56*EV_[8]*EV_[10]
                H_[9,7] = H_[7,9]
                H_[7,10] = V12*V34*V56*EV_[8]*EV_[9]
                H_[10,7] = H_[7,10]
                H_[8,9] = V12*V34*V56*EV_[7]*EV_[10]
                H_[9,8] = H_[8,9]
                H_[8,10] = V12*V34*V56*EV_[7]*EV_[9]
                H_[10,8] = H_[8,10]
                H_[9,10] = V12*V34*V56*V78
                H_[10,9] = H_[9,10]
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

