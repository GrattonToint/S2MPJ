function SENSORS(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : SENSORS
#    *********
# 
#    A problem arising from two-dimensional optimal sensor placement
# 
#    Source:
#    H. Zhang and X. Wang,
#    "Optimal sensor placement",
#    SIAM Review, vol. 35, p. 641, 1993.
# 
#    SIF input: Nick Gould, June 1994
# 
#    classification = "C-OUR2-AN-V-0"
# 
#    Number of unknowns
# 
#       Alternative values for the SIF file parameters:
# IE N                   2              $-PARAMETER
# IE N                   3              $-PARAMETER
# IE N                   10             $-PARAMETER
# IE N                   100            $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "SENSORS"

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
# IE N                   1000           $-PARAMETER
        v_["1"] = 1
        v_["RN"] = Float64(v_["N"])
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("THETA"*string(I),ix_)
            arrset(pb.xnames,iv,"THETA"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for J = Int64(v_["1"]):Int64(v_["N"])
            for I = Int64(v_["1"]):Int64(v_["N"])
                ig,ig_,_ = s2mpj_ii("S"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<>")
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
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        for I = Int64(v_["1"]):Int64(v_["N"])
            v_["RI"] = Float64(I)
            v_["I/N"] = v_["RI"]/v_["RN"]
            pb.x0[ix_["THETA"*string(I)]] = Float64(v_["I/N"])
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSINFUN", iet_)
        loaset(elftv,it,1,"THETAI")
        loaset(elftv,it,2,"THETAJ")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for J = Int64(v_["1"]):Int64(v_["N"])
            for I = Int64(v_["1"]):Int64(v_["N"])
                ename = "S"*string(I)*","*string(J)
                ie,ie_,newelt = s2mpj_ii(ename,ie_)
                if newelt > 0
                    arrset(pbm.elftype,ie,"eSINFUN")
                    arrset(ielftype,ie,iet_["eSINFUN"])
                end
                vname = "THETA"*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="THETAI",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "THETA"*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="THETAJ",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gmL2",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for ig in 1:ngrp
            arrset(pbm.grftype,ig,"gmL2")
        end
        for J = Int64(v_["1"]):Int64(v_["N"])
            for I = Int64(v_["1"]):Int64(v_["N"])
                ig = ig_["S"*string(I)*","*string(J)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["S"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel,1.)
            end
        end
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.A = spzeros(Float64,0,0)
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-OUR2-AN-V-0"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eSINFUN"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        TIMJ = EV_[1]-EV_[2]
        SI = sin(EV_[1])
        SJ = sin(EV_[2])
        SIMJ = sin(TIMJ)
        CI = cos(EV_[1])
        CJ = cos(EV_[2])
        CIMJ = cos(TIMJ)
        CJSIMJ = CJ*SIMJ-SJ*CIMJ
        CJCIMJ = CJ*CIMJ+SJ*SIMJ
        f_   = SI*SJ*SIMJ
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = SJ*(CI*SIMJ+SI*CIMJ)
            g_[2] = SI*CJSIMJ
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = 2.0*SJ*(CI*CIMJ-SI*SIMJ)
                H_[1,2] = CI*CJSIMJ+SI*CJCIMJ
                H_[2,1] = H_[1,2]
                H_[2,2] = -2.0*SI*CJCIMJ
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

    elseif action == "gmL2"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_= -GVAR_*GVAR_
        if nargout>1
            g_ = -GVAR_-GVAR_
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = -2.0
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

