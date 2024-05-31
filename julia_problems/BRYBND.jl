function BRYBND(action,args...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : BRYBND
#    *********
#    Broyden banded system of nonlinear equations, considered in the
#    least square sense.
# 
#    Source: problem 31 in
#    J.J. More', B.S. Garbow and K.E. Hillstrom,
#    "Testing Unconstrained Optimization Software",
#    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
# 
#    See also Buckley#73 (p. 41) and Toint#18
# 
#    SDIF input: Ph. Toint, Dec 1989.
# 
#    classification = "SUR2-AN-V-0"
# 
#    N is the number of equations and variables (variable).
# 
# IE N                   10             $-PARAMETER     original value
# IE N                   50             $-PARAMETER
# IE N                   100            $-PARAMETER
# IE N                   500            $-PARAMETER
# IE N                   1000           $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "BRYBND"

    if action == "setup"
        pbm          = PBM(name)
        pb           = PB(name)
        pb.sifpbname = "BRYBND"
        nargin       = length(args)
        pbm.call     = eval( Meta.parse( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        if nargin<1
            v_["N"] = 5000;  #  SIF file default value
        else
            v_["N"] = args[1];
        end
        if nargin<2
            v_["KAPPA1"] = 2.0;  #  SIF file default value
        else
            v_["KAPPA1"] = args[2];
        end
        if nargin<3
            v_["KAPPA2"] = 5.0;  #  SIF file default value
        else
            v_["KAPPA2"] = args[3];
        end
        if nargin<4
            v_["KAPPA3"] = 1.0;  #  SIF file default value
        else
            v_["KAPPA3"] = args[4];
        end
        if nargin<5
            v_["LB"] = 5;  #  SIF file default value
        else
            v_["LB"] = args[5];
        end
        if nargin<6
            v_["UB"] = 1;  #  SIF file default value
        else
            v_["UB"] = args[6];
        end
        v_["1"] = 1
        v_["MLB"] = -1*v_["LB"]
        v_["MUB"] = -1*v_["UB"]
        v_["LB+1"] = 1+v_["LB"]
        v_["N-UB"] = v_["N"]+v_["MUB"]
        v_["N-UB-1"] = -1+v_["N-UB"]
        v_["-KAPPA3"] = -1.0*v_["KAPPA3"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        xscale  = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2x_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["LB"])
            v_["I-1"] = -1+I
            v_["I+1"] = 1+I
            v_["I+UB"] = I+v_["UB"]
            for J = Int64(v_["1"]):Int64(v_["I-1"])
                ig,ig_,_ = s2x_ii("G"*string(I),ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["X"*string(J)]
                pbm.A[ig,iv] += v_["-KAPPA3"]
            end
            ig,ig_,_ = s2x_ii("G"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += v_["KAPPA1"]
            for J = Int64(v_["I+1"]):Int64(v_["I+UB"])
                ig,ig_,_ = s2x_ii("G"*string(I),ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["X"*string(J)]
                pbm.A[ig,iv] += v_["-KAPPA3"]
            end
        end
        for I = Int64(v_["LB+1"]):Int64(v_["N-UB-1"])
            v_["I-LB"] = I+v_["MLB"]
            v_["I-1"] = -1+I
            v_["I+1"] = 1+I
            v_["I+UB"] = I+v_["UB"]
            for J = Int64(v_["I-LB"]):Int64(v_["I-1"])
                ig,ig_,_ = s2x_ii("G"*string(I),ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["X"*string(J)]
                pbm.A[ig,iv] += v_["-KAPPA3"]
            end
            ig,ig_,_ = s2x_ii("G"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += v_["KAPPA1"]
            for J = Int64(v_["I+1"]):Int64(v_["I+UB"])
                ig,ig_,_ = s2x_ii("G"*string(I),ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["X"*string(J)]
                pbm.A[ig,iv] += v_["-KAPPA3"]
            end
        end
        for I = Int64(v_["N-UB"]):Int64(v_["N"])
            v_["I-LB"] = I+v_["MLB"]
            v_["I-1"] = -1+I
            v_["I+1"] = 1+I
            for J = Int64(v_["I-LB"]):Int64(v_["I-1"])
                ig,ig_,_ = s2x_ii("G"*string(I),ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["X"*string(J)]
                pbm.A[ig,iv] += v_["-KAPPA3"]
            end
            ig,ig_,_ = s2x_ii("G"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += v_["KAPPA1"]
            for J = Int64(v_["I+1"]):Int64(v_["N"])
                ig,ig_,_ = s2x_ii("G"*string(I),ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["X"*string(J)]
                pbm.A[ig,iv] += v_["-KAPPA3"]
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
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(1.0,pb.n,)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2x_ii( "SQ", iet_)
        loaset(elftv,it,1,"V")
        it,iet_,_ = s2x_ii( "CB", iet_)
        loaset(elftv,it,1,"V")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["N"])
            ename = "E"*string(I)
            ie,ie_,newelt = s2x_ii(ename,ie_)
            arrset(pbm.elftype,ie,"SQ")
            arrset(ielftype, ie, iet_["SQ"])
            vname = "X"*string(I)
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,1.0)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "Q"*string(I)
            ie,ie_,newelt = s2x_ii(ename,ie_)
            arrset(pbm.elftype,ie,"CB")
            arrset(ielftype, ie, iet_["CB"])
            vname = "X"*string(I)
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,1.0)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2x_ii("L2",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for ig in 1:ngrp
            arrset(pbm.grftype,ig,"L2")
        end
        for I = Int64(v_["1"]):Int64(v_["LB"])
            v_["I-1"] = -1+I
            v_["I+1"] = 1+I
            v_["I+UB"] = I+v_["UB"]
            for J = Int64(v_["1"]):Int64(v_["I-1"])
                ig = ig_["G"*string(I)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["E"*string(J)])
                loaset(pbm.grelw,ig,posel,v_["-KAPPA3"])
            end
            ig = ig_["G"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["Q"*string(I)])
            loaset(pbm.grelw,ig,posel,v_["KAPPA2"])
            for J = Int64(v_["I+1"]):Int64(v_["I+UB"])
                ig = ig_["G"*string(I)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["E"*string(J)])
                loaset(pbm.grelw,ig,posel,v_["-KAPPA3"])
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
                loaset(pbm.grelw,ig,posel,v_["-KAPPA3"])
            end
            ig = ig_["G"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E"*string(I)])
            loaset(pbm.grelw,ig,posel,v_["KAPPA2"])
            for J = Int64(v_["I+1"]):Int64(v_["I+UB"])
                ig = ig_["G"*string(I)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["E"*string(J)])
                loaset(pbm.grelw,ig,posel,v_["-KAPPA3"])
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
                loaset(pbm.grelw,ig,posel,v_["-KAPPA3"])
            end
            ig = ig_["G"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["Q"*string(I)])
            loaset(pbm.grelw,ig,posel,v_["KAPPA2"])
            for J = Int64(v_["I+1"]):Int64(v_["N"])
                ig = ig_["G"*string(I)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["E"*string(J)])
                loaset(pbm.grelw,ig,posel,v_["-KAPPA3"])
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTIONS %%%%%%%
        pb.pbclass = "SUR2-AN-V-0"
        return pb, pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "SQ"

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

    elseif action == "CB"

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

    elseif action == "L2"

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

    elseif action in  ["fx","fgx","fgHx","cx","cJx","cJHx","cIx","cIJx","cIJHx","cIJxv","fHxv","cJxv","Lxy","Lgxy","LgHxy","LIxy","LIgxy","LIgHxy","LHxyv","LIHxyv"]

        pbm = args[1]
        if pbm.name == name
            pbm.has_globs = [0,0]
            return s2x_eval(action,args...)
        else
            println("ERROR: please run "*name*" with action = setup")
            return ntuple(i->undef,args[end])
        end

    else
        println("ERROR: unknown action "*action*" requested from "*name*"%s.jl")
        return ntuple(i->undef,args[end])
    end

end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

