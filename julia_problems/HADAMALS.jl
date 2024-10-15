function HADAMALS(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HADAMALS
#    --------
# 
#    An attempt to find Hadamard matrices of order N.
# 
#    The problem is to find an N by N orthonormal matrix Q,
#    with column norms N, whose entries are plus or minus one.
# 
#    Source:  A suggestion by Alan Edelman (MIT).
# 
#    SIF input: Nick Gould, Nov 1993.
# 
#    classification = "C-OBR2-RN-V-V"
# 
#    The dimension of the matrix (=> N**2 variables).
# 
#       Alternative values for the SIF file parameters:
# IE N                   2              $-PARAMETER    original value
# IE N                   4              $-PARAMETER
# IE N                   6              $-PARAMETER
# IE N                   8              $-PARAMETER
# IE N                   10             $-PARAMETER
# IE N                   12             $-PARAMETER
# IE N                   14             $-PARAMETER
# IE N                   16             $-PARAMETER
# IE N                   18             $-PARAMETER
# IE N                   20             $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HADAMALS"

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
# IE N                   32             $-PARAMETER
# IE N                   64             $-PARAMETER
# IE N                   128            $-PARAMETER
# IE N                   256            $-PARAMETER
# IE N                   428            $-PARAMETER
        v_["1"] = 1
        v_["2"] = 2
        v_["RN"] = Float64(v_["N"])
        v_["N/2"] = trunc(Int,(v_["N"]/v_["2"]))
        v_["N/2+1"] = 1+v_["N/2"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for J = Int64(v_["1"]):Int64(v_["N"])
            for I = Int64(v_["1"]):Int64(v_["N"])
                iv,ix_,_ = s2mpj_ii("Q"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"Q"*string(I)*","*string(J))
            end
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for J = Int64(v_["1"]):Int64(v_["N"])
            for I = Int64(v_["1"]):Int64(J)
                ig,ig_,_ = s2mpj_ii("O"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<>")
            end
        end
        for J = Int64(v_["1"]):Int64(v_["N"])
            for I = Int64(v_["2"]):Int64(v_["N"])
                ig,ig_,_ = s2mpj_ii("S"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<>")
            end
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        for J = Int64(v_["1"]):Int64(v_["N"])
            pbm.gconst[ig_["O"*string(J)*","*string(J)]] = Float64(v_["RN"])
        end
        for J = Int64(v_["1"]):Int64(v_["N"])
            for I = Int64(v_["2"]):Int64(v_["N"])
                pbm.gconst[ig_["S"*string(I)*","*string(J)]] = Float64(1.0)
            end
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = fill(-1.0,pb.n)
        pb.xupper = fill(1.0,pb.n)
        for I = Int64(v_["1"]):Int64(v_["N/2"])
            pb.xlower[ix_["Q"*string(I)*","*string(Int64(v_["1"]))]] = 1.0
            pb.xupper[ix_["Q"*string(I)*","*string(Int64(v_["1"]))]] = 1.0
        end
        for I = Int64(v_["N/2+1"]):Int64(v_["N"])
            pb.xlower[ix_["Q"*string(I)*","*string(Int64(v_["1"]))]] = -1.0
            pb.xupper[ix_["Q"*string(I)*","*string(Int64(v_["1"]))]] = -1.0
        end
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        for J = Int64(v_["1"]):Int64(v_["N"])
            for I = Int64(v_["1"]):Int64(v_["N/2"])
                pb.x0[ix_["Q"*string(I)*","*string(J)]] = Float64(0.9)
            end
            for I = Int64(v_["N/2+1"]):Int64(v_["N"])
                pb.x0[ix_["Q"*string(I)*","*string(J)]] = Float64(-0.9)
            end
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSQR", iet_)
        loaset(elftv,it,1,"Q1")
        it,iet_,_ = s2mpj_ii( "en2PROD", iet_)
        loaset(elftv,it,1,"Q1")
        loaset(elftv,it,2,"Q2")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for J = Int64(v_["1"]):Int64(v_["N"])
            for I = Int64(v_["1"]):Int64(J)
                for K = Int64(v_["1"]):Int64(v_["N"])
                    ename = "O"*string(I)*","*string(J)*","*string(K)
                    ie,ie_,_  = s2mpj_ii(ename,ie_)
                    arrset(pbm.elftype,ie,"en2PROD")
                    arrset(ielftype,ie,iet_["en2PROD"])
                    vname = "Q"*string(K)*","*string(I)
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-1.0),Float64(1.0),nothing)
                    posev = findfirst(x->x=="Q1",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    vname = "Q"*string(K)*","*string(J)
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-1.0),Float64(1.0),nothing)
                    posev = findfirst(x->x=="Q2",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                end
            end
        end
        for J = Int64(v_["1"]):Int64(v_["N"])
            for I = Int64(v_["2"]):Int64(v_["N"])
                ename = "S"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eSQR")
                arrset(ielftype,ie,iet_["eSQR"])
                vname = "Q"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-1.0),Float64(1.0),nothing)
                posev = findfirst(x->x=="Q1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gL2",igt_)
        it,igt_,_ = s2mpj_ii("gLARGEL2",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for J = Int64(v_["1"]):Int64(v_["N"])
            for I = Int64(v_["1"]):Int64(J)
                ig = ig_["O"*string(I)*","*string(J)]
                arrset(pbm.grftype,ig,"gL2")
                for K = Int64(v_["1"]):Int64(v_["N"])
                    ig = ig_["O"*string(I)*","*string(J)]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["O"*string(I)*","*string(J)*","*string(K)])
                    loaset(pbm.grelw,ig,posel,1.)
                end
            end
        end
        for J = Int64(v_["1"]):Int64(v_["N"])
            for I = Int64(v_["2"]):Int64(v_["N"])
                ig = ig_["S"*string(I)*","*string(J)]
                arrset(pbm.grftype,ig,"gLARGEL2")
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["S"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel,1.)
            end
        end
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.A = spzeros(Float64,0,0)
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-OBR2-RN-V-V"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm


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
                H_[1,1] = 2.0e+0
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "en2PROD"

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
                H_[1,2] = 1.0e+0
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
    elseif action == "g_globs"

        pbm = args[1]
        arrset(pbm.gfpar,1,1.0e+0)    # this is  FACTOR
        return pbm

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
                H_ = 2.0e+0
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "gLARGEL2"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_= pbm.gfpar[1]*GVAR_*GVAR_
        if nargout>1
            g_ = 2.0e+0*pbm.gfpar[1]*GVAR_
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = 2.0e+0*pbm.gfpar[1]
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
            pbm.has_globs = [0,1]
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

