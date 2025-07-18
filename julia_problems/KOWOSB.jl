function KOWOSB(action::String,args::Union{Any}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : KOWOSB
#    *********
# 
#    A problem arising in the analysis of kinetic data for an enzyme
#    reaction, known under the name of Kowalik and Osborne problem
#    in 4 variables.
# 
#    Source:  Problem 15 in
#    J.J. More', B.S. Garbow and K.E. Hillstrom,
#    "Testing Unconstrained Optimization Software",
#    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
# 
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "C-CSUR2-MN-4-0"
# 
#    This function  is a nonlinear least squares with 11 groups.  Each
#    group has a linear and a nonlinear element.
# 
#    Number of groups
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 21 VI 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "KOWOSB"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling KOWOSB.")
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
        v_["M"] = 11
        v_["1"] = 1
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
        iv,ix_,_ = s2mpj_ii("X4",ix_)
        arrset(pb.xnames,iv,"X4")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype = String[]
        for I = Int64(v_["1"]):Int64(v_["M"])
            ig,ig_,_ = s2mpj_ii("G"*string(I),ig_)
            arrset(gtype,ig,"<>")
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        pbm.gconst[ig_["G1"]] = Float64(0.1957)
        pbm.gconst[ig_["G2"]] = Float64(0.1947)
        pbm.gconst[ig_["G3"]] = Float64(0.1735)
        pbm.gconst[ig_["G4"]] = Float64(0.1600)
        pbm.gconst[ig_["G5"]] = Float64(0.0844)
        pbm.gconst[ig_["G6"]] = Float64(0.0627)
        pbm.gconst[ig_["G7"]] = Float64(0.0456)
        pbm.gconst[ig_["G8"]] = Float64(0.0342)
        pbm.gconst[ig_["G9"]] = Float64(0.0323)
        pbm.gconst[ig_["G10"]] = Float64(0.0235)
        pbm.gconst[ig_["G11"]] = Float64(0.0246)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.x0[ix_["X1"]] = Float64(0.25)
        pb.x0[ix_["X2"]] = Float64(0.39)
        pb.x0[ix_["X3"]] = Float64(0.415)
        pb.x0[ix_["X4"]] = Float64(0.39)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eKWO", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        loaset(elftv,it,3,"V3")
        loaset(elftv,it,4,"V4")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"U")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["M"])
            ename = "E"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eKWO")
            arrset(ielftype,ie,iet_["eKWO"])
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
        end
        ename = "E1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        posep = findfirst(x->x=="U",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(4.0))
        ename = "E2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        posep = findfirst(x->x=="U",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(2.0))
        ename = "E3"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        posep = findfirst(x->x=="U",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        ename = "E4"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        posep = findfirst(x->x=="U",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.5))
        ename = "E5"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        posep = findfirst(x->x=="U",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.25))
        ename = "E6"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        posep = findfirst(x->x=="U",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.167))
        ename = "E7"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        posep = findfirst(x->x=="U",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.125))
        ename = "E8"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        posep = findfirst(x->x=="U",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.1))
        ename = "E9"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        posep = findfirst(x->x=="U",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.0833))
        ename = "E10"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        posep = findfirst(x->x=="U",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.0714))
        ename = "E11"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        posep = findfirst(x->x=="U",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.0624))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gL2",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["M"])
            ig = ig_["G"*string(I)]
            arrset(pbm.grftype,ig,"gL2")
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E"*string(I)])
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        pb.objlower = 0.0
#    Solution
# LO SOLTN               0.00102734
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-CSUR2-MN-4-0"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eKWO"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        USQ = pbm.elpar[iel_][1]*pbm.elpar[iel_][1]
        B1 = USQ+pbm.elpar[iel_][1]*EV_[2]
        B2 = USQ+pbm.elpar[iel_][1]*EV_[3]+EV_[4]
        B2SQ = B2*B2
        B2CB = B2*B2SQ
        UV1 = pbm.elpar[iel_][1]*EV_[1]
        UB1 = pbm.elpar[iel_][1]*B1
        T1 = B1/B2SQ
        T2 = 2.0/B2CB
        f_   = EV_[1]*B1/B2
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = B1/B2
            g_[2] = UV1/B2
            g_[3] = -UV1*T1
            g_[4] = -EV_[1]*T1
            if nargout>2
                H_ = zeros(Float64,4,4)
                H_[1,2] = pbm.elpar[iel_][1]/B2
                H_[2,1] = H_[1,2]
                H_[1,3] = -UB1/B2SQ
                H_[3,1] = H_[1,3]
                H_[1,4] = -T1
                H_[4,1] = H_[1,4]
                H_[2,3] = -UV1*pbm.elpar[iel_][1]/B2SQ
                H_[3,2] = H_[2,3]
                H_[2,4] = -UV1/B2SQ
                H_[4,2] = H_[2,4]
                H_[3,3] = T2*UV1*UB1
                H_[3,4] = T2*UV1*B1
                H_[4,3] = H_[3,4]
                H_[4,4] = T2*EV_[1]*B1
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

