function PT(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : PT
#    *********
# 
#    A nonlinear programming formulation of a discretization of
#    a nonlinear minimax problem.
# 
#    The problem is
# 
#        min  max phi(x,w), for all w in the interval I.
#         x    w
# 
#    I is discretized, and the problem solved over the
#    discrete points.
# 
#    Nonlinear programming formulation
#        min   u     s.t.  u - phi >= 0
#        x,u
# 
#    Specific problem: I = [0,1]
#    phi(x,w) = (2w^2-1)x + w(1-w)(1-x)
# 
#    Source: E. Polak and A. L. Tits,
#    "A recursive quadratic programming algorithm for semi-infinite
#     optimization problems",
#    Appl. Math. Optim. 8, 1982, pp 325-349.
# 
#    SIF input: Nick Gould, February, 1994.
# 
#    classification = "C-LLR2-AN-2-V"
# 
#    Discretization
# 
# IE M                   2
# IE M                   100
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "PT"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["M"] = 500
        v_["LOWER"] = 0.0
        v_["UPPER"] = 1.0
        v_["ONE"] = 1.0
        v_["0"] = 0
        v_["DIFF"] = v_["UPPER"]-v_["LOWER"]
        v_["RM"] = Float64(v_["M"])
        v_["H"] = v_["DIFF"]/v_["RM"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        iv,ix_,_ = s2mpj_ii("U",ix_)
        arrset(pb.xnames,iv,"U")
        iv,ix_,_ = s2mpj_ii("X",ix_)
        arrset(pb.xnames,iv,"X")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["U"]
        pbm.A[ig,iv] += Float64(1.0)
        for I = Int64(v_["0"]):Int64(v_["M"])
            v_["RI"] = Float64(I)
            v_["W"] = v_["RI"]*v_["H"]
            v_["W"] = v_["W"]+v_["LOWER"]
            v_["1-W"] = 1.0-v_["W"]
            v_["W(1-W)"] = v_["W"]*v_["1-W"]
            v_["W**2"] = v_["W"]*v_["W"]
            v_["2W**2"] = 2.0*v_["W**2"]
            v_["2W**2-1"] = v_["2W**2"]-v_["ONE"]
            v_["XCOEFF"] = v_["W(1-W)"]-v_["2W**2-1"]
            ig,ig_,_ = s2mpj_ii("LO"*string(I),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"LO"*string(I))
            iv = ix_["U"]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["X"]
            pbm.A[ig,iv] += Float64(v_["XCOEFF"])
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
        for I = Int64(v_["0"]):Int64(v_["M"])
            v_["RI"] = Float64(I)
            v_["W"] = v_["RI"]*v_["H"]
            v_["W"] = v_["W"]+v_["LOWER"]
            v_["1-W"] = 1.0-v_["W"]
            v_["W-W**2"] = v_["W"]*v_["1-W"]
            pbm.gconst[ig_["LO"*string(I)]] = Float64(v_["W-W**2"])
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+pb.neq+1:pb.m] = zeros(Float64,pb.nge)
        pb.cupper[1:pb.nge] = fill(Inf,pb.nge)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = collect(1:length(pbm.congrps))
        pb.pbclass = "C-LLR2-AN-2-V"
        pb.x0          = zeros(Float64,pb.n)
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm


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

