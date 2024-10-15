function BROYDNBDLS(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : BROYDNBDLS
#    *********
#    Broyden banded system of nonlinear equations, considered in the
#    least square sense.
# 
#    Source: problem 31 in
#    J.J. More', B.S. Garbow and K.E. Hillstrom,
#    "Testing Unconstrained Optimization Software",
#    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
# 
#    See also Buckley#73 and Toint#18
#    SIF input: Ph. Toint, Dec 1989.
#    Least-squares version: Nick Gould, Oct 2015
# 
#    classification = "C-SUR2-AN-V-0"
# 
#    N is the number of equations and variables (variable).
# 
#       Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER     original value
# IE N                   50             $-PARAMETER
# IE N                   100            $-PARAMETER
# IE N                   500            $-PARAMETER
# IE N                   1000           $-PARAMETER
# IE N                   5000           $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "BROYDNBDLS"

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
# IE N                   10000          $-PARAMETER
        if nargin<2
            v_["KAPPA1"] = Float64(2.0);  #  SIF file default value
        else
            v_["KAPPA1"] = Float64(args[2]);
        end
        if nargin<3
            v_["KAPPA2"] = Float64(5.0);  #  SIF file default value
        else
            v_["KAPPA2"] = Float64(args[3]);
        end
        if nargin<4
            v_["KAPPA3"] = Float64(1.0);  #  SIF file default value
        else
            v_["KAPPA3"] = Float64(args[4]);
        end
        if nargin<5
            v_["LB"] = Int64(5);  #  SIF file default value
        else
            v_["LB"] = Int64(args[5]);
        end
        if nargin<6
            v_["UB"] = Int64(1);  #  SIF file default value
        else
            v_["UB"] = Int64(args[6]);
        end
        v_["1"] = 1
        v_["MLB"] = -1*v_["LB"]
        v_["MUB"] = -1*v_["UB"]
        v_["LB+1"] = 1+v_["LB"]
        v_["N-UB"] = v_["N"]+v_["MUB"]
        v_["N-UB-1"] = -1+v_["N-UB"]
        v_["-KAPPA3"] = -1.0*v_["KAPPA3"]
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
        for I = Int64(v_["1"]):Int64(v_["LB"])
            v_["I-1"] = -1+I
            v_["I+1"] = 1+I
            v_["I+UB"] = I+v_["UB"]
            for J = Int64(v_["1"]):Int64(v_["I-1"])
                ig,ig_,_ = s2mpj_ii("G"*string(I),ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["X"*string(J)]
                pbm.A[ig,iv] += Float64(v_["-KAPPA3"])
            end
            ig,ig_,_ = s2mpj_ii("G"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(v_["KAPPA1"])
            for J = Int64(v_["I+1"]):Int64(v_["I+UB"])
                ig,ig_,_ = s2mpj_ii("G"*string(I),ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["X"*string(J)]
                pbm.A[ig,iv] += Float64(v_["-KAPPA3"])
            end
        end
        for I = Int64(v_["LB+1"]):Int64(v_["N-UB-1"])
            v_["I-LB"] = I+v_["MLB"]
            v_["I-1"] = -1+I
            v_["I+1"] = 1+I
            v_["I+UB"] = I+v_["UB"]
            for J = Int64(v_["I-LB"]):Int64(v_["I-1"])
                ig,ig_,_ = s2mpj_ii("G"*string(I),ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["X"*string(J)]
                pbm.A[ig,iv] += Float64(v_["-KAPPA3"])
            end
            ig,ig_,_ = s2mpj_ii("G"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(v_["KAPPA1"])
            for J = Int64(v_["I+1"]):Int64(v_["I+UB"])
                ig,ig_,_ = s2mpj_ii("G"*string(I),ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["X"*string(J)]
                pbm.A[ig,iv] += Float64(v_["-KAPPA3"])
            end
        end
        for I = Int64(v_["N-UB"]):Int64(v_["N"])
            v_["I-LB"] = I+v_["MLB"]
            v_["I-1"] = -1+I
            v_["I+1"] = 1+I
            for J = Int64(v_["I-LB"]):Int64(v_["I-1"])
                ig,ig_,_ = s2mpj_ii("G"*string(I),ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["X"*string(J)]
                pbm.A[ig,iv] += Float64(v_["-KAPPA3"])
            end
            ig,ig_,_ = s2mpj_ii("G"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(v_["KAPPA1"])
            for J = Int64(v_["I+1"]):Int64(v_["N"])
                ig,ig_,_ = s2mpj_ii("G"*string(I),ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["X"*string(J)]
                pbm.A[ig,iv] += Float64(v_["-KAPPA3"])
            end
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(1.0),pb.n)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSQ", iet_)
        loaset(elftv,it,1,"V")
        it,iet_,_ = s2mpj_ii( "eCB", iet_)
        loaset(elftv,it,1,"V")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["N"])
            ename = "E"*string(I)
            ie,ie_,newelt = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQ")
            arrset(ielftype,ie,iet_["eSQ"])
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "Q"*string(I)
            ie,ie_,newelt = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eCB")
            arrset(ielftype,ie,iet_["eCB"])
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
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
        for I = Int64(v_["1"]):Int64(v_["LB"])
            v_["I-1"] = -1+I
            v_["I+1"] = 1+I
            v_["I+UB"] = I+v_["UB"]
            for J = Int64(v_["1"]):Int64(v_["I-1"])
                ig = ig_["G"*string(I)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["E"*string(J)])
                loaset(pbm.grelw,ig,posel,Float64(v_["-KAPPA3"]))
            end
            ig = ig_["G"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["Q"*string(I)])
            loaset(pbm.grelw,ig,posel,Float64(v_["KAPPA2"]))
            for J = Int64(v_["I+1"]):Int64(v_["I+UB"])
                ig = ig_["G"*string(I)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["E"*string(J)])
                loaset(pbm.grelw,ig,posel,Float64(v_["-KAPPA3"]))
            end
        end
        for I = Int64(v_["LB+1"]):Int64(v_["N-UB-1"])
            v_["I-LB"] = I+v_["MLB"]
            v_["I-1"] = -1+I
            v_["I+1"] = 1+I
            v_["I+UB"] = I+v_["UB"]
            for J = Int64(v_["I-LB"]):Int64(v_["I-1"])
                ig = ig_["G"*string(I)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["Q"*string(J)])
                loaset(pbm.grelw,ig,posel,Float64(v_["-KAPPA3"]))
            end
            ig = ig_["G"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E"*string(I)])
            loaset(pbm.grelw,ig,posel,Float64(v_["KAPPA2"]))
            for J = Int64(v_["I+1"]):Int64(v_["I+UB"])
                ig = ig_["G"*string(I)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["E"*string(J)])
                loaset(pbm.grelw,ig,posel,Float64(v_["-KAPPA3"]))
            end
        end
        for I = Int64(v_["N-UB"]):Int64(v_["N"])
            v_["I-LB"] = I+v_["MLB"]
            v_["I-1"] = -1+I
            v_["I+1"] = 1+I
            for J = Int64(v_["I-LB"]):Int64(v_["I-1"])
                ig = ig_["G"*string(I)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["E"*string(J)])
                loaset(pbm.grelw,ig,posel,Float64(v_["-KAPPA3"]))
            end
            ig = ig_["G"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["Q"*string(I)])
            loaset(pbm.grelw,ig,posel,Float64(v_["KAPPA2"]))
            for J = Int64(v_["I+1"]):Int64(v_["N"])
                ig = ig_["G"*string(I)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["E"*string(J)])
                loaset(pbm.grelw,ig,posel,Float64(v_["-KAPPA3"]))
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
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

    elseif action == "eSQ"

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

    elseif action == "eCB"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[1]*EV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 3.0*EV_[1]*EV_[1]
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 6.0*EV_[1]
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

