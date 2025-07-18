function SPECAN(action::String,args::Union{Any}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : SPECAN
#    *********
# 
#    Source: a problem in spectral analysis suggested
#    by J. Eriksson and P. Lindstrom in "A Parallel Algorithm
#    for Bound Constrained Nonlinear Least Squares", UMEA TR S-901 87
# 
#    SIF input: Michael Ferris, July 1993
# 
#    classification = "C-CSBR2-AN-V-0"
# 
#    Number of Gaussians
# 
#       Alternative values for the SIF file parameters:
# IE K                   1              $-PARAMETER
# IE K                   2              $-PARAMETER
# IE K                   3              $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 21 VI 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "SPECAN"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling SPECAN.")
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
        if nargin<1
            v_["K"] = Int64(3);  #  SIF file default value
        else
            v_["K"] = Int64(args[1]);
        end
        v_["N"] = 3
        v_["M"] = 5000
        v_["RealM"] = Float64(v_["M"])
        v_["H"] = 25.0/v_["RealM"]
        v_["1"] = 1
        v_["2"] = 2
        v_["3"] = 3
        v_["ONE"] = 1.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        irA   = Int64[]
        icA   = Int64[]
        valA  = Float64[]
        for p = Int64(v_["1"]):Int64(v_["K"])
            for j = Int64(v_["1"]):Int64(v_["N"])
                iv,ix_,_ = s2mpj_ii("X"*string(p)*","*string(j),ix_)
                arrset(pb.xnames,iv,"X"*string(p)*","*string(j))
            end
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype = String[]
        for p = Int64(v_["1"]):Int64(v_["K"])
            for I = Int64(v_["1"]):Int64(v_["M"])
                ig,ig_,_ = s2mpj_ii("OBJ"*string(p)*","*string(I),ig_)
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
        v_["SOLN1,1"] = 19.0
        v_["SOLN1,2"] = 4.2
        v_["SOLN1,3"] = 1.2
        v_["SOLN2,1"] = 8.0
        v_["SOLN2,2"] = 2.5
        v_["SOLN2,3"] = 4.6
        v_["SOLN3,1"] = 10.0
        v_["SOLN3,2"] = 2.0
        v_["SOLN3,3"] = 2.6
        for I = Int64(v_["1"]):Int64(v_["M"])
            v_["RI"] = Float64(I)
            v_["IH"] = v_["H"]*v_["RI"]
            v_["TI"] = v_["ONE"]+v_["IH"]
            v_["Differ"] = v_["TI"]-v_["SOLN1,2"]
            v_["Numer"] = v_["Differ"]*v_["Differ"]
            v_["Denom"] = v_["SOLN1,3"]*v_["SOLN1,3"]
            v_["Differ"] = v_["Numer"]/v_["Denom"]
            v_["Ratio"] = 0.0-v_["Differ"]
            v_["ERat"] = exp(v_["Ratio"])
            v_["Yi1"] = v_["SOLN1,1"]*v_["ERat"]
            v_["Differ"] = v_["TI"]-v_["SOLN2,2"]
            v_["Numer"] = v_["Differ"]*v_["Differ"]
            v_["Denom"] = v_["SOLN2,3"]*v_["SOLN2,3"]
            v_["Differ"] = v_["Numer"]/v_["Denom"]
            v_["Ratio"] = 0.0-v_["Differ"]
            v_["ERat"] = exp(v_["Ratio"])
            v_["Yi2"] = v_["SOLN2,1"]*v_["ERat"]
            v_["Differ"] = v_["TI"]-v_["SOLN3,2"]
            v_["Numer"] = v_["Differ"]*v_["Differ"]
            v_["Denom"] = v_["SOLN3,3"]*v_["SOLN3,3"]
            v_["Differ"] = v_["Numer"]/v_["Denom"]
            v_["Ratio"] = 0.0-v_["Differ"]
            v_["ERat"] = exp(v_["Ratio"])
            v_["Yi3"] = v_["SOLN3,1"]*v_["ERat"]
            for p = Int64(v_["1"]):Int64(v_["K"])
                pbm.gconst[ig_["OBJ"*string(p)*","*string(I)]] = Float64(v_["Yi"*string(p)])
            end
        end
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        v_["LOWER1,1"] = 15.0
        v_["LOWER1,2"] = 3.5
        v_["LOWER1,3"] = 0.3
        v_["LOWER2,1"] = 5.0
        v_["LOWER2,2"] = 2.2
        v_["LOWER2,3"] = 2.6
        v_["LOWER3,1"] = 5.0
        v_["LOWER3,2"] = 1.2
        v_["LOWER3,3"] = 1.3
        v_["UPPER1,1"] = 31.0
        v_["UPPER1,2"] = 6.3
        v_["UPPER1,3"] = 3.7
        v_["UPPER2,1"] = 15.0
        v_["UPPER2,2"] = 5.3
        v_["UPPER2,3"] = 6.2
        v_["UPPER3,1"] = 14.0
        v_["UPPER3,2"] = 3.3
        v_["UPPER3,3"] = 2.8
        for p = Int64(v_["1"]):Int64(v_["K"])
            for j = Int64(v_["1"]):Int64(v_["N"])
                pb.xlower[ix_["X"*string(p)*","*string(j)]]  = (
                      v_["LOWER"*string(p)*","*string(j)])
                pb.xupper[ix_["X"*string(p)*","*string(j)]]  = (
                      v_["UPPER"*string(p)*","*string(j)])
            end
        end
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        v_["START1,1"] = 25.0
        v_["START1,2"] = 5.2
        v_["START1,3"] = 3.2
        v_["START2,1"] = 7.0
        v_["START2,2"] = 4.1
        v_["START2,3"] = 3.6
        v_["START3,1"] = 11.6
        v_["START3,2"] = 1.9
        v_["START3,3"] = 2.2
        for p = Int64(v_["1"]):Int64(v_["K"])
            for j = Int64(v_["1"]):Int64(v_["N"])
                pb.x0[ix_["X"*string(p)*","*string(j)]]  = (
                      Float64(v_["START"*string(p)*","*string(j)]))
            end
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eEXPSQ", iet_)
        loaset(elftv,it,1,"U")
        loaset(elftv,it,2,"V")
        loaset(elftv,it,3,"W")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"T")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for p = Int64(v_["1"]):Int64(v_["K"])
            for I = Int64(v_["1"]):Int64(v_["M"])
                ename = "E"*string(p)*","*string(I)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eEXPSQ")
                arrset(ielftype,ie,iet_["eEXPSQ"])
                vname = "X"*string(p)*","*string(Int64(v_["1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="U",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(p)*","*string(Int64(v_["2"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(p)*","*string(Int64(v_["3"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="W",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                v_["RI"] = Float64(I)
                v_["IH"] = v_["H"]*v_["RI"]
                v_["TI"] = v_["ONE"]+v_["IH"]
                posep = findfirst(x->x=="T",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["TI"]))
            end
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gSQUARE",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for ig in 1:ngrp
            arrset(pbm.grftype,ig,"gSQUARE")
        end
        for p = Int64(v_["1"]):Int64(v_["K"])
            for I = Int64(v_["1"]):Int64(v_["M"])
                ig = ig_["OBJ"*string(p)*","*string(I)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["E"*string(p)*","*string(I)])
                loaset(pbm.grelw,ig,posel,1.)
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-CSBR2-AN-V-0"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eEXPSQ"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        R = (pbm.elpar[iel_][1]-EV_[2])^2
        S = EV_[3]^2
        E = exp(-R/S)
        f_   = EV_[1]*E
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = E
            g_[2] = 2.0*(pbm.elpar[iel_][1]-EV_[2])*EV_[1]*E/S
            g_[3] = 2.0*R*EV_[1]*E/(S*EV_[3])
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,2] = 2.0*(pbm.elpar[iel_][1]-EV_[2])*E/S
                H_[2,1] = H_[1,2]
                H_[1,3] = 2.0*R*E/(S*EV_[3])
                H_[3,1] = H_[1,3]
                H_[2,2] = (2.0*EV_[1]*E/S)*(2.0*R/S-1.0)
                H_[2,3] = 4.0*(pbm.elpar[iel_][1]-EV_[2])*EV_[1]*E/(S*EV_[3])*(R/S-1.0)
                H_[3,2] = H_[2,3]
                H_[3,3] = 2.0*R*EV_[1]*E/(S^3)*(2.0*R-3.0*S)
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

    elseif action == "gSQUARE"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_= 5.0e-1*GVAR_*GVAR_
        if nargout>1
            g_ = GVAR_
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = 1.0e+0
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

