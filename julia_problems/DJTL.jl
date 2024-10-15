function DJTL(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DJTL
#    *********
# 
#    Source: modified version of problem 19 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
#    that is meant to simulate the Lagrangian barrier objective function
#    for particular values of the shifts and multipliers
# 
#    SIF input: A.R. Conn August 1993
# 
#    classification = "C-OUR2-AN-2-0"
# 
#    Define multipliers and shifts
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "DJTL"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["LL1"] = 1.0
        v_["LL2"] = 1.0
        v_["LL3"] = 1.0
        v_["LL4"] = 1.0
        v_["SL1"] = 1.0
        v_["SL2"] = 1.0
        v_["SL3"] = 1.0
        v_["SL4"] = 1.0
        v_["LU1"] = 1.0
        v_["LU2"] = 1.0
        v_["LU3"] = 1.0
        v_["LU4"] = 1.0
        v_["SU1"] = 1.0
        v_["SU2"] = 1.0
        v_["SU3"] = 1.0
        v_["SU4"] = 1.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        iv,ix_,_ = s2mpj_ii("X1",ix_)
        arrset(pb.xnames,iv,"X1")
        iv,ix_,_ = s2mpj_ii("X2",ix_)
        arrset(pb.xnames,iv,"X2")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("CONU1",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("CONL1",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("CONU2",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("CONL2",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("BNDU1",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("BNDL1",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("BNDU2",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("BNDL2",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(1.0)
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        pbm.gconst[ig_["CONU1"]] = Float64(-200.0)
        pbm.gconst[ig_["CONL1"]] = Float64(100.0)
        pbm.gconst[ig_["CONU2"]] = Float64(0.0)
        pbm.gconst[ig_["CONL2"]] = Float64(-82.81)
        pbm.gconst[ig_["BNDU1"]] = Float64(-100.0)
        pbm.gconst[ig_["BNDL1"]] = Float64(13.0)
        pbm.gconst[ig_["BNDU2"]] = Float64(-100.0)
        pbm.gconst[ig_["BNDL2"]] = Float64(0.0)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.x0[ix_["X1"]] = Float64(15.0)
        pb.x0[ix_["X2"]] = Float64(6.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eCBm10", iet_)
        loaset(elftv,it,1,"V1")
        it,iet_,_ = s2mpj_ii( "eCBm20", iet_)
        loaset(elftv,it,1,"V1")
        it,iet_,_ = s2mpj_ii( "eSQm5", iet_)
        loaset(elftv,it,1,"V1")
        it,iet_,_ = s2mpj_ii( "eSQm6", iet_)
        loaset(elftv,it,1,"V1")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "E1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eCBm10")
        arrset(ielftype,ie,iet_["eCBm10"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eCBm20")
        arrset(ielftype,ie,iet_["eCBm20"])
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E3"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQm5")
        arrset(ielftype,ie,iet_["eSQm5"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E4"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQm5")
        arrset(ielftype,ie,iet_["eSQm5"])
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E5"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQm6")
        arrset(ielftype,ie,iet_["eSQm6"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gLOG",igt_)
        it,igt_,_ = s2mpj_ii("gLOG",igt_)
        grftp = Vector{Vector{String}}()
        loaset(grftp,it,1,"P1")
        it,igt_,_ = s2mpj_ii("gLOG",igt_)
        loaset(grftp,it,2,"P2")
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["OBJ"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E1"])
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E2"])
        loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_["CONL1"]
        arrset(pbm.grftype,ig,"gLOG")
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E3"])
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E4"])
        loaset(pbm.grelw,ig,posel, 1.)
        posgp = findfirst(x->x=="P1",grftp[igt_[pbm.grftype[ig]]])
        loaset(pbm.grpar,ig,posgp,Float64(v_["SL1"]))
        posgp = findfirst(x->x=="P2",grftp[igt_[pbm.grftype[ig]]])
        loaset(pbm.grpar,ig,posgp,Float64(v_["LL1"]))
        ig = ig_["CONU1"]
        arrset(pbm.grftype,ig,"gLOG")
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E3"])
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E4"])
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        posgp = findfirst(x->x=="P1",grftp[igt_[pbm.grftype[ig]]])
        loaset(pbm.grpar,ig,posgp,Float64(v_["SU1"]))
        posgp = findfirst(x->x=="P2",grftp[igt_[pbm.grftype[ig]]])
        loaset(pbm.grpar,ig,posgp,Float64(v_["LU1"]))
        ig = ig_["CONL2"]
        arrset(pbm.grftype,ig,"gLOG")
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E4"])
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E5"])
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        posgp = findfirst(x->x=="P1",grftp[igt_[pbm.grftype[ig]]])
        loaset(pbm.grpar,ig,posgp,Float64(v_["SL2"]))
        posgp = findfirst(x->x=="P2",grftp[igt_[pbm.grftype[ig]]])
        loaset(pbm.grpar,ig,posgp,Float64(v_["LL2"]))
        ig = ig_["CONU2"]
        arrset(pbm.grftype,ig,"gLOG")
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E4"])
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E5"])
        loaset(pbm.grelw,ig,posel, 1.)
        posgp = findfirst(x->x=="P1",grftp[igt_[pbm.grftype[ig]]])
        loaset(pbm.grpar,ig,posgp,Float64(v_["SU2"]))
        posgp = findfirst(x->x=="P2",grftp[igt_[pbm.grftype[ig]]])
        loaset(pbm.grpar,ig,posgp,Float64(v_["LU2"]))
        ig = ig_["BNDL1"]
        arrset(pbm.grftype,ig,"gLOG")
        posgp = findfirst(x->x=="P1",grftp[igt_[pbm.grftype[ig]]])
        loaset(pbm.grpar,ig,posgp,Float64(v_["SL3"]))
        posgp = findfirst(x->x=="P2",grftp[igt_[pbm.grftype[ig]]])
        loaset(pbm.grpar,ig,posgp,Float64(v_["LL3"]))
        ig = ig_["BNDU1"]
        arrset(pbm.grftype,ig,"gLOG")
        posgp = findfirst(x->x=="P1",grftp[igt_[pbm.grftype[ig]]])
        loaset(pbm.grpar,ig,posgp,Float64(v_["SU3"]))
        posgp = findfirst(x->x=="P2",grftp[igt_[pbm.grftype[ig]]])
        loaset(pbm.grpar,ig,posgp,Float64(v_["LU3"]))
        ig = ig_["BNDL2"]
        arrset(pbm.grftype,ig,"gLOG")
        posgp = findfirst(x->x=="P1",grftp[igt_[pbm.grftype[ig]]])
        loaset(pbm.grpar,ig,posgp,Float64(v_["SL4"]))
        posgp = findfirst(x->x=="P2",grftp[igt_[pbm.grftype[ig]]])
        loaset(pbm.grpar,ig,posgp,Float64(v_["LL4"]))
        ig = ig_["BNDU2"]
        arrset(pbm.grftype,ig,"gLOG")
        posgp = findfirst(x->x=="P1",grftp[igt_[pbm.grftype[ig]]])
        loaset(pbm.grpar,ig,posgp,Float64(v_["SU4"]))
        posgp = findfirst(x->x=="P2",grftp[igt_[pbm.grftype[ig]]])
        loaset(pbm.grpar,ig,posgp,Float64(v_["LU4"]))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               -8951.54472
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-OUR2-AN-2-0"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eCBm10"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        DIF = EV_[1]-10.0
        f_   = DIF^3
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 3.0*DIF*DIF
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 6.0*DIF
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eCBm20"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        DIF = EV_[1]-20.0
        f_   = DIF^3
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 3.0*DIF*DIF
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 6.0*DIF
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eSQm5"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        DIF = EV_[1]-5.0
        f_   = DIF^2
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 2.0*DIF
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

    elseif action == "eSQm6"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        DIF = EV_[1]-6.0
        f_   = DIF^2
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 2.0*DIF
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

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    elseif action == "gLOG"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        APP1 = GVAR_+pbm.grpar[igr_][1]
        P1P2 = pbm.grpar[igr_][1]*pbm.grpar[igr_][2]
        ARG0 = APP1<=0.0
        BIG = 1.0000e+10
        if ARG0
            FF = BIG*GVAR_^2
        end
        if !ARG0
            FF = -P1P2*log(APP1)
        end
        if ARG0
            GG = 2.0*BIG*GVAR_
        end
        if !ARG0
            GG = -P1P2/APP1
        end
        if ARG0
            HH = 2.0*BIG
        end
        if !ARG0
            HH = P1P2/APP1^2
        end
        f_= FF
        if nargout>1
            g_ = GG
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = HH
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

