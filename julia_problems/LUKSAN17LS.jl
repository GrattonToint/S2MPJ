function LUKSAN17LS(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : LUKSAN17LS
#    *********
# 
#    Problem 17 (sparse trigonometric) in the paper
# 
#      L. Luksan
#      Hybrid methods in large sparse nonlinear least squares
#      J. Optimization Theory & Applications 89(3) 575-595 (1996)
# 
#    SIF input: Nick Gould, June 2017.
# 
#    least-squares version
# 
#    classification = "C-SUR2-AN-V-0"
# 
#   seed for dimensions
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "LUKSAN17LS"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["S"] = 49
        v_["N"] = 2*v_["S"]
        v_["N"] = 2+v_["N"]
        v_["M"] = 4*v_["S"]
        v_["1"] = 1
        v_["2"] = 2
        v_["3"] = 3
        v_["4"] = 4
        v_["ONE"] = 1.0
        v_["Y1"] = 30.6
        v_["Y2"] = 72.2
        v_["Y3"] = 124.4
        v_["Y4"] = 187.4
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
        for I = Int64(v_["1"]):Int64(v_["M"])
            ig,ig_,_ = s2mpj_ii("E"*string(I),ig_)
            arrset(gtype,ig,"<>")
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        v_["K"] = 1
        for J = Int64(v_["1"]):Int64(v_["S"])
            pbm.gconst[ig_["E"*string(Int64(v_["K"]))]] = Float64(v_["Y1"])
            v_["K"] = 1+v_["K"]
            pbm.gconst[ig_["E"*string(Int64(v_["K"]))]] = Float64(v_["Y2"])
            v_["K"] = 1+v_["K"]
            pbm.gconst[ig_["E"*string(Int64(v_["K"]))]] = Float64(v_["Y3"])
            v_["K"] = 1+v_["K"]
            pbm.gconst[ig_["E"*string(Int64(v_["K"]))]] = Float64(v_["Y4"])
            v_["K"] = 1+v_["K"]
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        for I = Int64(v_["1"]):Int64(v_["4"]):Int64(v_["N"])
            pb.x0[ix_["X"*string(I)]] = Float64(-0.8)
        end
        for I = Int64(v_["2"]):Int64(v_["4"]):Int64(v_["N"])
            pb.x0[ix_["X"*string(I)]] = Float64(1.2)
        end
        for I = Int64(v_["3"]):Int64(v_["4"]):Int64(v_["N"])
            pb.x0[ix_["X"*string(I)]] = Float64(-1.2)
        end
        for I = Int64(v_["4"]):Int64(v_["4"]):Int64(v_["N"])
            pb.x0[ix_["X"*string(I)]] = Float64(0.8)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eACOSX", iet_)
        loaset(elftv,it,1,"X")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"A")
        it,iet_,_ = s2mpj_ii( "eASINX", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftp,it,1,"A")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for Q = Int64(v_["1"]):Int64(v_["4"])
            v_["RQ"] = Float64(Q)
            v_["RQ2"] = v_["RQ"]*v_["RQ"]
            v_["K"] = 1
            v_["I"] = 0
            for J = Int64(v_["1"]):Int64(v_["S"])
                v_["I+Q"] = v_["I"]+Q
                for L = Int64(v_["1"]):Int64(v_["4"])
                    v_["RL"] = Float64(L)
                    v_["RL2"] = v_["RL"]*v_["RL"]
                    v_["A"] = v_["RL"]*v_["RQ2"]
                    v_["A"] = -1.0*v_["A"]
                    ename = "S"*string(Int64(v_["K"]))*","*string(Q)
                    ie,ie_,_  = s2mpj_ii(ename,ie_)
                    arrset(pbm.elftype,ie,"eASINX")
                    arrset(ielftype,ie,iet_["eASINX"])
                    ename = "S"*string(Int64(v_["K"]))*","*string(Q)
                    ie,ie_,_  = s2mpj_ii(ename,ie_)
                    vname = "X"*string(Int64(v_["I+Q"]))
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="X",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    ename = "S"*string(Int64(v_["K"]))*","*string(Q)
                    ie,ie_,_  = s2mpj_ii(ename,ie_)
                    posep = findfirst(x->x=="A",elftp[ielftype[ie]])
                    loaset(pbm.elpar,ie,posep,Float64(v_["A"]))
                    v_["A"] = v_["RL2"]*v_["RQ"]
                    ename = "C"*string(Int64(v_["K"]))*","*string(Q)
                    ie,ie_,_  = s2mpj_ii(ename,ie_)
                    arrset(pbm.elftype,ie,"eACOSX")
                    arrset(ielftype,ie,iet_["eACOSX"])
                    ename = "C"*string(Int64(v_["K"]))*","*string(Q)
                    ie,ie_,_  = s2mpj_ii(ename,ie_)
                    vname = "X"*string(Int64(v_["I+Q"]))
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="X",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    ename = "C"*string(Int64(v_["K"]))*","*string(Q)
                    ie,ie_,_  = s2mpj_ii(ename,ie_)
                    posep = findfirst(x->x=="A",elftp[ielftype[ie]])
                    loaset(pbm.elpar,ie,posep,Float64(v_["A"]))
                    v_["K"] = 1+v_["K"]
                end
                v_["I"] = 2+v_["I"]
            end
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
        for K = Int64(v_["1"]):Int64(v_["M"])
            for Q = Int64(v_["1"]):Int64(v_["4"])
                ig = ig_["E"*string(K)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["S"*string(K)*","*string(Q)])
                loaset(pbm.grelw,ig,posel,1.)
                posel = posel+1
                loaset(pbm.grelt,ig,posel,ie_["C"*string(K)*","*string(Q)])
                loaset(pbm.grelw,ig,posel, 1.)
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN                0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.A = spzeros(Float64,0,0)
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

    elseif action == "eASINX"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        ASINX = pbm.elpar[iel_][1]*sin(EV_[1])
        f_   = ASINX
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = pbm.elpar[iel_][1]*cos(EV_[1])
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = -ASINX
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eACOSX"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        ACOSX = pbm.elpar[iel_][1]*cos(EV_[1])
        f_   = ACOSX
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = -pbm.elpar[iel_][1]*sin(EV_[1])
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = -ACOSX
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

