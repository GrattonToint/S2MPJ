function EXPQUAD(action::String,args::Union{Any}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : EXPQUAD
#    *********
# 
#    A problem with mixed exponential and quadratic terms.
# 
#    SIF input: Ph. Toint, 1992.
#               minor correction by Ph. Shott, Jan 1995.
# 
#    classification = "C-COBR2-AN-V-V"
# 
#       Alternative values for the SIF file parameters:
# IE N                   12             $-PARAMETER
# IE M                   6              $-PARAMETER
# 
# IE N                   120            $-PARAMETER     original value
# IE M                   10             $-PARAMETER     original value
# 
# IE N                   1200           $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 22 VII 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "EXPQUAD"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling EXPQUAD.")
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
            v_["N"] = Int64(12);  #  SIF file default value
        else
            v_["N"] = Int64(args[1]);
        end
# IE M                   100            $-PARAMETER
        if nargin<2
            v_["M"] = Int64(6);  #  SIF file default value
        else
            v_["M"] = Int64(args[2]);
        end
        v_["0"] = 0
        v_["1"] = 1
        v_["2"] = 2
        v_["M+1"] = 1+v_["M"]
        v_["N-1"] = -1+v_["N"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        irA   = Int64[]
        icA   = Int64[]
        valA  = Float64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype = String[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            v_["RI"] = Float64(I)
            v_["C"] = - 10.0*v_["RI"]
            ig,ig_,_ = s2mpj_ii("OBJ",ig_)
            arrset(gtype,ig,"<>")
            push!(irA,ig)
            push!(icA,ix_["X"*string(I)])
            push!(valA,Float64(v_["C"]))
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        for I = Int64(v_["1"]):Int64(v_["M"])
            pb.xlower[ix_["X"*string(I)]] = 0.0
            pb.xupper[ix_["X"*string(I)]] = 10.0
        end
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eEXP", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"P")
        it,iet_,_ = s2mpj_ii( "eQUAD", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        v_["RM"] = Float64(v_["M"])
        for I = Int64(v_["1"]):Int64(v_["M"])
            v_["RI"] = Float64(I)
            v_["C"] = v_["RI"]/v_["RM"]
            v_["I+1"] = I+v_["1"]
            ename = "E"*string(I)
            ie,ie_,newelt = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eEXP")
            arrset(ielftype,ie,iet_["eEXP"])
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["I+1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="P",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["C"]))
        end
        for I = Int64(v_["M+1"]):Int64(v_["N-1"])
            ename = "E"*string(I)
            ie,ie_,newelt = s2mpj_ii(ename,ie_)
            if newelt > 0
                arrset(pbm.elftype,ie,"eQUAD")
                arrset(ielftype,ie,iet_["eQUAD"])
            end
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["N"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            ig = ig_["OBJ"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E"*string(I)])
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-COBR2-AN-V-V"
        pb.x0          = zeros(Float64,pb.n)
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eEXP"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        F = exp(0.1*pbm.elpar[iel_][1]*EV_[1]*EV_[2])
        f_   = F
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 0.1*pbm.elpar[iel_][1]*EV_[2]*F
            g_[2] = 0.1*pbm.elpar[iel_][1]*EV_[1]*F
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = (0.1*pbm.elpar[iel_][1]*EV_[2])^2*F
                H_[2,2] = (0.1*pbm.elpar[iel_][1]*EV_[1])^2*F
                H_[2,1] = (0.1+0.01*pbm.elpar[iel_][1]*EV_[1]*EV_[2])*F*pbm.elpar[iel_][1]
                H_[1,2] = H_[2,1]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eQUAD"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = 4.0*EV_[1]*EV_[1]+2.0*EV_[2]*EV_[2]+EV_[1]*EV_[2]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 8.0*EV_[1]+EV_[2]
            g_[2] = 4.0*EV_[2]+EV_[1]
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = 8.0
                H_[1,2] = 1.0
                H_[2,1] = H_[1,2]
                H_[2,2] = 4.0
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

