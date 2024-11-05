function LUKSAN11(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : LUKSAN11
#    *********
# 
#    Problem 11 (chained serpentine) in the paper
# 
#      L. Luksan
#      Hybrid methods in large sparse nonlinear least squares
#      J. Optimization Theory & Applications 89(3) 575-595 (1996)
# 
#    SIF input: Nick Gould, June 2017.
# 
#    classification = "C-CNOR2-AN-V-V"
# 
#   seed for dimensions
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 17 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "LUKSAN11"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["S"] = 99
        v_["N"] = v_["S"]
        v_["N"] = 1+v_["N"]
        v_["M"] = 2*v_["S"]
        v_["1"] = 1
        v_["2"] = 2
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
        v_["K"] = 1
        for I = Int64(v_["1"]):Int64(v_["S"])
            v_["I+1"] = 1+I
            ig,ig_,_ = s2mpj_ii("E"*string(Int64(v_["K"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"E"*string(Int64(v_["K"])))
            iv = ix_["X"*string(Int64(v_["I+1"]))]
            pbm.A[ig,iv] += Float64(-10.0)
            v_["K"] = 1+v_["K"]
            ig,ig_,_ = s2mpj_ii("E"*string(Int64(v_["K"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"E"*string(Int64(v_["K"])))
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            v_["K"] = 1+v_["K"]
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        legrps = findall(x->x=="<=",gtype)
        eqgrps = findall(x->x=="==",gtype)
        gegrps = findall(x->x==">=",gtype)
        pb.nle = length(legrps)
        pb.neq = length(eqgrps)
        pb.nge = length(gegrps)
        pb.m   = pb.nle+pb.neq+pb.nge
        pbm.congrps = [[legrps;eqgrps];gegrps]
        pb.nob = ngrp-pb.m
        pbm.objgrps = findall(x->x=="<>",gtype)
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        for I = Int64(v_["1"]):Int64(v_["2"]):Int64(v_["M"])
            v_["I+1"] = 1+I
            pbm.gconst[ig_["E"*string(Int64(v_["I+1"]))]] = Float64(1.0)
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        for I = Int64(v_["1"]):Int64(v_["N"])
            pb.x0[ix_["X"*string(I)]] = Float64(-0.8)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eL1", iet_)
        loaset(elftv,it,1,"V")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        v_["K"] = 1
        for I = Int64(v_["1"]):Int64(v_["S"])
            ename = "L"*string(Int64(v_["K"]))
            ie,ie_,newelt = s2mpj_ii(ename,ie_)
            if newelt > 0
                arrset(pbm.elftype,ie,"eL1")
                arrset(ielftype,ie,iet_["eL1"])
            end
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            v_["K"] = 2+v_["K"]
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        v_["K"] = 1
        for I = Int64(v_["1"]):Int64(v_["S"])
            ig = ig_["E"*string(Int64(v_["K"]))]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["L"*string(Int64(v_["K"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            v_["K"] = 2+v_["K"]
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN               0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-CNOR2-AN-V-V"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eL1"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        D = 1.0e0+EV_[1]*EV_[1]
        f_   = 2.0e1*EV_[1]/D
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 2.0e1*(1.0e0-EV_[1]*EV_[1])/D^2
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = -4.0e1*EV_[1]*(3.0e0-EV_[1]*EV_[1])/D^3
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

