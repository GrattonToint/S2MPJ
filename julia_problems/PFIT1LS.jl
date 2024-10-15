function PFIT1LS(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : PFIT1LS
#    *********
# 
#    The problem is to fit a model containing a pole, given data
#    for values, first and second derivatives at two distinct points.
#    This is a least-squares version of problem PFIT1.
# 
#    The problem is not convex.
# 
#    SIF input: Ph. Toint, March 1994.
#               Lower bound on H added, Nov 2002.
# 
#    classification = "C-SBR2-AN-3-0"
# 
#    Problem data
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "PFIT1LS"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["CF"] = -8.0
        v_["CG"] = -18.66666666
        v_["CH"] = -23.11111111
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        iv,ix_,_ = s2mpj_ii("A",ix_)
        arrset(pb.xnames,iv,"A")
        iv,ix_,_ = s2mpj_ii("R",ix_)
        arrset(pb.xnames,iv,"R")
        iv,ix_,_ = s2mpj_ii("H",ix_)
        arrset(pb.xnames,iv,"H")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("EF",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("EG",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("EH",ig_)
        arrset(gtype,ig,"<>")
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        pbm.gconst[ig_["EF"]] = Float64(v_["CF"])
        pbm.gconst[ig_["EG"]] = Float64(v_["CG"])
        pbm.gconst[ig_["EH"]] = Float64(v_["CH"])
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        pb.xlower[ix_["H"]] = -0.5
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.x0[ix_["A"]] = Float64(1.0)
        pb.x0[ix_["R"]] = Float64(0.0)
        pb.x0[ix_["H"]] = Float64(1.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eT1", iet_)
        loaset(elftv,it,1,"AA")
        loaset(elftv,it,2,"RR")
        loaset(elftv,it,3,"XX")
        it,iet_,_ = s2mpj_ii( "eT2", iet_)
        loaset(elftv,it,1,"AA")
        loaset(elftv,it,2,"RR")
        loaset(elftv,it,3,"XX")
        it,iet_,_ = s2mpj_ii( "eT3", iet_)
        loaset(elftv,it,1,"AA")
        loaset(elftv,it,2,"RR")
        loaset(elftv,it,3,"XX")
        it,iet_,_ = s2mpj_ii( "eT4", iet_)
        loaset(elftv,it,1,"AA")
        loaset(elftv,it,2,"RR")
        loaset(elftv,it,3,"XX")
        it,iet_,_ = s2mpj_ii( "eT5", iet_)
        loaset(elftv,it,1,"AA")
        loaset(elftv,it,2,"RR")
        loaset(elftv,it,3,"XX")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "EA"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eT3")
        arrset(ielftype,ie,iet_["eT3"])
        vname = "A"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="AA",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "R"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RR",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "H"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EB"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eT2")
        arrset(ielftype,ie,iet_["eT2"])
        vname = "A"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="AA",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "R"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RR",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "H"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EC"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eT1")
        arrset(ielftype,ie,iet_["eT1"])
        vname = "A"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="AA",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "R"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RR",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "H"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "ED"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eT4")
        arrset(ielftype,ie,iet_["eT4"])
        vname = "A"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="AA",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "R"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RR",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "H"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EE"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eT5")
        arrset(ielftype,ie,iet_["eT5"])
        vname = "A"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="AA",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "R"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RR",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "H"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
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
        ig = ig_["EF"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EA"])
        loaset(pbm.grelw,ig,posel,Float64(-0.5))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["EC"])
        loaset(pbm.grelw,ig,posel,Float64(1.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["ED"])
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        ig = ig_["EG"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EA"])
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["EB"])
        loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_["EH"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EE"])
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution at ( 1.0, 3.0 , 2.0 )
# LO SOLTN               0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.A = spzeros(Float64,0,0)
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-SBR2-AN-3-0"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eT1"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[2]*EV_[3]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[2]*EV_[3]
            g_[2] = EV_[1]*EV_[3]
            g_[3] = EV_[1]*EV_[2]
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,2] = EV_[3]
                H_[2,1] = H_[1,2]
                H_[1,3] = EV_[2]
                H_[3,1] = H_[1,3]
                H_[2,3] = EV_[1]
                H_[3,2] = H_[2,3]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eT2"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        A1 = EV_[1]+1.0
        Y = 1.0+EV_[3]
        LOGY = log(Y)
        C = Y^(-A1)
        CC = C/Y
        CCC = CC/Y
        B = 1.0-C
        BA = LOGY*C
        BX = A1*CC
        BAA = -LOGY*LOGY*C
        BAX = -LOGY*BX+CC
        BXX = -A1*(A1+1.0)*CCC
        ARX = EV_[1]*EV_[2]*EV_[3]
        f_   = ARX*B
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[2]*EV_[3]*B+ARX*BA
            g_[2] = EV_[1]*EV_[3]*B
            g_[3] = EV_[1]*EV_[2]*B+ARX*BX
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,1] = 2.0*EV_[2]*EV_[3]*BA+ARX*BAA
                H_[1,2] = EV_[3]*B+EV_[1]*EV_[3]*BA
                H_[2,1] = H_[1,2]
                H_[1,3] = EV_[2]*B+EV_[2]*EV_[3]*BX+EV_[1]*EV_[2]*BA+ARX*BAX
                H_[3,1] = H_[1,3]
                H_[2,3] = EV_[1]*B+EV_[1]*EV_[3]*BX
                H_[3,2] = H_[2,3]
                H_[3,3] = 2.0*EV_[1]*EV_[2]*BX+ARX*BXX
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eT3"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*(EV_[1]+1.0)*EV_[2]*EV_[3]*EV_[3]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = (2.0*EV_[1]+1.0)*EV_[2]*EV_[3]*EV_[3]
            g_[2] = EV_[1]*(EV_[1]+1.0)*EV_[3]*EV_[3]
            g_[3] = 2.0*EV_[1]*(EV_[1]+1.0)*EV_[2]*EV_[3]
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,1] = 2.0*EV_[2]*EV_[3]*EV_[3]
                H_[1,2] = (2.0*EV_[1]+1.0)*EV_[3]*EV_[3]
                H_[2,1] = H_[1,2]
                H_[1,3] = 2.0*(2.0*EV_[1]+1.0)*EV_[2]*EV_[3]
                H_[3,1] = H_[1,3]
                H_[2,3] = 2.0*EV_[1]*(EV_[1]+1.0)*EV_[3]
                H_[3,2] = H_[2,3]
                H_[3,3] = 2.0*EV_[1]*(EV_[1]+1.0)*EV_[2]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eT4"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        Y = 1.0+EV_[3]
        LOGY = log(Y)
        C = Y^(-EV_[1])
        CC = C/Y
        CCC = CC/Y
        B = 1.0-C
        BA = LOGY*C
        BX = EV_[1]*CC
        BAA = -LOGY*LOGY*C
        BAX = -LOGY*BX+CC
        BXX = -EV_[1]*(EV_[1]+1.0)*CCC
        f_   = EV_[2]*B
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[2]*BA
            g_[2] = B
            g_[3] = EV_[2]*BX
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,1] = EV_[2]*BAA
                H_[1,2] = BA
                H_[2,1] = H_[1,2]
                H_[1,3] = EV_[2]*BAX
                H_[3,1] = H_[1,3]
                H_[2,3] = BX
                H_[3,2] = H_[2,3]
                H_[3,3] = EV_[2]*BXX
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eT5"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        A1 = EV_[1]+2.0
        Y = 1.0+EV_[3]
        LOGY = log(Y)
        C = Y^(-A1)
        CC = C/Y
        CCC = CC/Y
        B = 1.0-C
        BA = LOGY*C
        BX = A1*CC
        BAA = -LOGY*LOGY*C
        BAX = -LOGY*BX+CC
        BXX = -A1*(A1+1.0)*CCC
        D = EV_[1]*(EV_[1]+1.0)*EV_[2]*EV_[3]*EV_[3]
        DA = (2.0*EV_[1]+1.0)*EV_[2]*EV_[3]*EV_[3]
        DR = EV_[1]*(EV_[1]+1.0)*EV_[3]*EV_[3]
        DX = 2.0*EV_[1]*(EV_[1]+1.0)*EV_[2]*EV_[3]
        DAA = 2.0*EV_[2]*EV_[3]*EV_[3]
        DAR = (2.0*EV_[1]+1.0)*EV_[3]*EV_[3]
        DAX = 2.0*(2.0*EV_[1]+1.0)*EV_[2]*EV_[3]
        DRX = 2.0*EV_[1]*(EV_[1]+1.0)*EV_[3]
        DXX = 2.0*EV_[1]*(EV_[1]+1.0)*EV_[2]
        f_   = D*B
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = DA*B+D*BA
            g_[2] = DR*B
            g_[3] = DX*B+D*BX
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,1] = DAA*B+2.0*DA*BA+D*BAA
                H_[1,2] = DAR*B+DR*BA
                H_[2,1] = H_[1,2]
                H_[1,3] = DAX*B+DA*BX+DX*BA+D*BAX
                H_[3,1] = H_[1,3]
                H_[2,3] = DRX*B+DR*BX
                H_[3,2] = H_[2,3]
                H_[3,3] = DXX*B+2.0*DX*BX+D*BXX
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

