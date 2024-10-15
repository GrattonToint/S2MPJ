function LEVYMONT7(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : LEVYMONT7
#    *********
#    A global optimization example due to Levy & Montalvo 
#    This problem is one of the parameterised set LEVYMONT5-LEVYMONT10
# 
#    Source:  problem 7 in
# 
#    A. V. Levy and A. Montalvo
#    "The Tunneling Algorithm for the Global Minimization of Functions"
#    SIAM J. Sci. Stat. Comp. 6(1) 1985 15:29 
#    https://doi.org/10.1137/0906002
# 
#    SIF input: Nick Gould, August 2021
# 
#    classification = "C-SBR2-AY-4-0"
# 
#    N is the number of variables
# 
#       Alternative values for the SIF file parameters:
# IE N                   4              $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "LEVYMONT7"

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
            v_["N"] = Int64(5);  #  SIF file default value
        else
            v_["N"] = Int64(args[1]);
        end
        v_["A"] = 1.0
        v_["K"] = 10.0
        v_["L"] = 0.25
        v_["C"] = 0.75
        v_["1"] = 1
        v_["2"] = 2
        v_["PI/4"] = atan(1.0)
        v_["PI"] = 4.0*v_["PI/4"]
        v_["RN"] = Float64(v_["N"])
        v_["A-C"] = v_["A"]-v_["C"]
        v_["PI/N"] = v_["PI"]/v_["RN"]
        v_["KPI/N"] = v_["K"]*v_["PI/N"]
        v_["ROOTKPI/N"] = sqrt(v_["KPI/N"])
        v_["N/PI"] = v_["RN"]/v_["PI"]
        v_["N/KPI"] = v_["N/PI"]/v_["K"]
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
        for I = Int64(v_["1"]):Int64(v_["N"])
            ig,ig_,_ = s2mpj_ii("Q"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(v_["L"])
            arrset(pbm.gscale,ig,Float64(v_["N/PI"]))
            ig,ig_,_ = s2mpj_ii("N"*string(I),ig_)
            arrset(gtype,ig,"<>")
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        for I = Int64(v_["1"]):Int64(v_["N"])
            pbm.gconst[ig_["Q"*string(I)]] = Float64(v_["A-C"])
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = fill(-10.0,pb.n)
        pb.xupper = fill(10.0,pb.n)
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(8.0),pb.n)
        pb.x0[ix_["X1"]] = Float64(-8.0)
        pb.x0[ix_["X2"]] = Float64(8.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eS2", iet_)
        loaset(elftv,it,1,"X")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"L")
        loaset(elftp,it,2,"C")
        it,iet_,_ = s2mpj_ii( "ePS2", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Z")
        loaset(elftp,it,1,"L")
        loaset(elftp,it,2,"C")
        loaset(elftp,it,3,"A")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "E"*string(Int64(v_["1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eS2")
        arrset(ielftype,ie,iet_["eS2"])
        ename = "E"*string(Int64(v_["1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["1"]))
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-10.0),Float64(10.0),Float64(8.0)))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        posep = findfirst(x->x=="L",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["L"]))
        ename = "E"*string(Int64(v_["1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        posep = findfirst(x->x=="C",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["C"]))
        for I = Int64(v_["2"]):Int64(v_["N"])
            v_["I-1"] = I-v_["1"]
            ename = "E"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ePS2")
            arrset(ielftype,ie,iet_["ePS2"])
            vname = "X"*string(I)
            iv,ix_,pb  = (
                  s2mpj_nlx(vname,ix_,pb,1,Float64(-10.0),Float64(10.0),Float64(8.0)))
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["I-1"]))
            iv,ix_,pb  = (
                  s2mpj_nlx(vname,ix_,pb,1,Float64(-10.0),Float64(10.0),Float64(8.0)))
            posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="L",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["L"]))
            posep = findfirst(x->x=="C",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["C"]))
            posep = findfirst(x->x=="A",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["A"]))
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gL2",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            ig = ig_["Q"*string(I)]
            arrset(pbm.grftype,ig,"gL2")
            ig = ig_["N"*string(I)]
            arrset(pbm.grftype,ig,"gL2")
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E"*string(I)])
            loaset(pbm.grelw,ig,posel,Float64(v_["ROOTKPI/N"]))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN               0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-SBR2-AY-4-0"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "e_globs"

        pbm = args[1]
        arrset(pbm.efpar,1,4.0*atan(1.0e0))
        return pbm

    elseif action == "eS2"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        PIL = pbm.efpar[1]*pbm.elpar[iel_][1]
        V = PIL*EV_[1]+pbm.efpar[1]*pbm.elpar[iel_][2]
        SINV = sin(V)
        COSV = cos(V)
        f_   = SINV
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = PIL*COSV
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = -PIL*PIL*SINV
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "ePS2"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        PIL = pbm.efpar[1]*pbm.elpar[iel_][1]
        U = pbm.elpar[iel_][1]*EV_[2]+pbm.elpar[iel_][2]-pbm.elpar[iel_][3]
        V = PIL*EV_[1]+pbm.efpar[1]*pbm.elpar[iel_][2]
        SINV = sin(V)
        COSV = cos(V)
        f_   = U*SINV
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = PIL*U*COSV
            g_[2] = pbm.elpar[iel_][1]*SINV
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = -PIL*PIL*U*SINV
                H_[1,2] = pbm.elpar[iel_][1]*PIL*COSV
                H_[2,1] = H_[1,2]
                H_[2,2] = 0.0
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
            g_ = 2.0*GVAR_
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
            pbm.has_globs = [1,0]
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

