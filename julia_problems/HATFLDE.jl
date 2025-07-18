function HATFLDE(action::String,args::Union{Any}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HATFLDE
#    *********
# 
#    An exponential fitting  test problem from the OPTIMA user manual.
# 
#    Source:
#    "The OPTIMA user manual (issue No.8, p. 37)",
#    Numerical Optimization Centre, Hatfield Polytechnic (UK), 1989.
# 
#    SIF input: Ph. Toint, May 1990.
# 
#    classification = "C-CSUR2-AN-3-0"
# 
#    Problem data
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 21 VI 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HATFLDE"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling HATFLDE.")
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
        v_["NG"] = 21
        v_["T1"] = 0.3
        v_["T2"] = 0.35
        v_["T3"] = 0.4
        v_["T4"] = 0.45
        v_["T5"] = 0.5
        v_["T6"] = 0.55
        v_["T7"] = 0.6
        v_["T8"] = 0.65
        v_["T9"] = 0.7
        v_["T10"] = 0.75
        v_["T11"] = 0.8
        v_["T12"] = 0.85
        v_["T13"] = 0.9
        v_["T14"] = 0.95
        v_["T15"] = 1.0
        v_["T16"] = 1.05
        v_["T17"] = 1.1
        v_["T18"] = 1.15
        v_["T19"] = 1.2
        v_["T20"] = 1.25
        v_["T21"] = 1.3
        v_["Z1"] = 1.561
        v_["Z2"] = 1.473
        v_["Z3"] = 1.391
        v_["Z4"] = 1.313
        v_["Z5"] = 1.239
        v_["Z6"] = 1.169
        v_["Z7"] = 1.103
        v_["Z8"] = 1.04
        v_["Z9"] = 0.981
        v_["Z10"] = 0.925
        v_["Z11"] = 0.8721
        v_["Z12"] = 0.8221
        v_["Z13"] = 0.7748
        v_["Z14"] = 0.73
        v_["Z15"] = 0.6877
        v_["Z16"] = 0.6477
        v_["Z17"] = 0.6099
        v_["Z18"] = 0.5741
        v_["Z19"] = 0.5403
        v_["Z20"] = 0.5084
        v_["Z21"] = 0.4782
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
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype = String[]
        for I = Int64(v_["1"]):Int64(v_["NG"])
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
        for I = Int64(v_["1"]):Int64(v_["NG"])
            v_["-ZI"] = -1.0*v_["Z"*string(I)]
            pbm.gconst[ig_["G"*string(I)]] = Float64(v_["-ZI"])
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.x0[ix_["X1"]] = Float64(1.0)
        pb.x0[ix_["X2"]] = Float64(-1.0)
        pb.x0[ix_["X3"]] = Float64(0.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eXPEXP", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"T")
        it,iet_,_ = s2mpj_ii( "ePEXP", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftp,it,1,"T")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["NG"])
            ename = "A"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eXPEXP")
            arrset(ielftype,ie,iet_["eXPEXP"])
            vname = "X1"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="T",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["T"*string(I)]))
            ename = "B"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ePEXP")
            arrset(ielftype,ie,iet_["ePEXP"])
            vname = "X3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="T",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["T"*string(I)]))
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
        for I = Int64(v_["1"]):Int64(v_["NG"])
            ig = ig_["G"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["A"*string(I)])
            loaset(pbm.grelw,ig,posel,Float64(-1.0))
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["B"*string(I)])
            loaset(pbm.grelw,ig,posel, 1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN               5.120377D-07
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-CSUR2-AN-3-0"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "ePEXP"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        EX = exp(pbm.elpar[iel_][1]*EV_[1])
        f_   = EX
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = pbm.elpar[iel_][1]*EX
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = pbm.elpar[iel_][1]*pbm.elpar[iel_][1]*EX
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eXPEXP"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        EX = exp(pbm.elpar[iel_][1]*EV_[2])
        f_   = EV_[1]*EX
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EX
            g_[2] = pbm.elpar[iel_][1]*EV_[1]*EX
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = pbm.elpar[iel_][1]*EX
                H_[2,1] = H_[1,2]
                H_[2,2] = pbm.elpar[iel_][1]*pbm.elpar[iel_][1]*EV_[1]*EX
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

