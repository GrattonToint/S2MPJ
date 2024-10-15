function OPTCNTRL(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : OPTCNTRL
#    *********
#    An optimal control problem
# 
#    Source:
#    B. Murtagh and M. Saunders,
#    Mathematical Programming studies 16, pp 84-117,
#    (example 5.11)
# 
#    SIF input: Nick Gould, June 1990.
# 
#    classification = "C-QQR2-AN-32-20"
# 
#    useful parameters
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "OPTCNTRL"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["T"] = 10
        v_["T-1"] = -1+v_["T"]
        v_["0"] = 0
        v_["1"] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for t = Int64(v_["0"]):Int64(v_["T"])
            iv,ix_,_ = s2mpj_ii("x"*string(t),ix_)
            arrset(pb.xnames,iv,"x"*string(t))
            iv,ix_,_ = s2mpj_ii("y"*string(t),ix_)
            arrset(pb.xnames,iv,"y"*string(t))
        end
        for t = Int64(v_["0"]):Int64(v_["T-1"])
            iv,ix_,_ = s2mpj_ii("u"*string(t),ix_)
            arrset(pb.xnames,iv,"u"*string(t))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        for t = Int64(v_["0"]):Int64(v_["T-1"])
            v_["t+1"] = 1+t
            ig,ig_,_ = s2mpj_ii("B"*string(t),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"B"*string(t))
            iv = ix_["x"*string(Int64(v_["t+1"]))]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["x"*string(t)]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["y"*string(t)]
            pbm.A[ig,iv] += Float64(-0.2)
            ig,ig_,_ = s2mpj_ii("C"*string(t),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"C"*string(t))
            iv = ix_["y"*string(Int64(v_["t+1"]))]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["y"*string(t)]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["x"*string(t)]
            pbm.A[ig,iv] += Float64(0.004)
            iv = ix_["u"*string(t)]
            pbm.A[ig,iv] += Float64(-0.2)
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
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        for t = Int64(v_["0"]):Int64(v_["T-1"])
            pb.xlower[ix_["x"*string(t)]] = -Inf
            pb.xupper[ix_["x"*string(t)]] = +Inf
            pb.xlower[ix_["y"*string(t)]] = -1.0
            pb.xlower[ix_["u"*string(t)]] = -0.2
            pb.xupper[ix_["u"*string(t)]] = 0.2
        end
        pb.xlower[ix_["x"*string(Int64(v_["0"]))]] = 10.0
        pb.xupper[ix_["x"*string(Int64(v_["0"]))]] = 10.0
        pb.xlower[ix_["y"*string(Int64(v_["0"]))]] = 0.0
        pb.xupper[ix_["y"*string(Int64(v_["0"]))]] = 0.0
        pb.xlower[ix_["y"*string(Int64(v_["T"]))]] = 0.0
        pb.xupper[ix_["y"*string(Int64(v_["T"]))]] = 0.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        for t = Int64(v_["1"]):Int64(v_["T-1"])
            if haskey(ix_,"y"*string(t))
                pb.x0[ix_["y"*string(t)]] = Float64(-1.0)
            else
                pb.y0[findfirst(x->x==ig_["y"*string(t)],pbm.congrps)] = Float64(-1.0)
            end
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSQR", iet_)
        loaset(elftv,it,1,"X")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for t = Int64(v_["0"]):Int64(v_["T"])
            ename = "o"*string(t)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
            vname = "x"*string(t)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        for t = Int64(v_["0"]):Int64(v_["T-1"])
            ename = "c"*string(t)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
            vname = "y"*string(t)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for t = Int64(v_["0"]):Int64(v_["T"])
            ig = ig_["OBJ"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["o"*string(t)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(0.5))
        end
        for t = Int64(v_["0"]):Int64(v_["T-1"])
            ig = ig_["C"*string(t)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["c"*string(t)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(0.01))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        pb.objlower = 0.0
#    Solution
# LO SOLTN               549.9999869
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
        pb.pbclass = "C-QQR2-AN-32-20"
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

    elseif action == "eSQR"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[1]+EV_[1]
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

