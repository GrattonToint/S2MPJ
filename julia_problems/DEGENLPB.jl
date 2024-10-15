function DEGENLPB(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DEGENLPB
#    *********
# 
#    A linear program with some degeneracy.
# 
#    Source:
#    T.C.T. Kotiah and D.I. Steinberg,
#    "Occurences of cycling and other phenomena arising in a class of
#    linear programming models",
#    Communications of the ACM, vol. 20, pp. 107-112, 1977.
# 
#    SIF input: Ph. Toint, Aug 1990.
# 
#    classification = "C-LLR2-AN-20-15"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "DEGENLPB"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["N"] = 20
        v_["M"] = 15
        v_["1"] = 1
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
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(-0.01)
        iv = ix_["X3"]
        pbm.A[ig,iv] += Float64(-33.333)
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(-100.0)
        iv = ix_["X5"]
        pbm.A[ig,iv] += Float64(-0.01)
        iv = ix_["X6"]
        pbm.A[ig,iv] += Float64(-33.343)
        iv = ix_["X7"]
        pbm.A[ig,iv] += Float64(-100.01)
        iv = ix_["X8"]
        pbm.A[ig,iv] += Float64(-33.333)
        iv = ix_["X9"]
        pbm.A[ig,iv] += Float64(-133.33)
        iv = ix_["X10"]
        pbm.A[ig,iv] += Float64(-100.0)
        ig,ig_,_ = s2mpj_ii("C1",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C1")
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(2.0)
        iv = ix_["X3"]
        pbm.A[ig,iv] += Float64(2.0)
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(2.0)
        iv = ix_["X5"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X6"]
        pbm.A[ig,iv] += Float64(2.0)
        iv = ix_["X7"]
        pbm.A[ig,iv] += Float64(2.0)
        iv = ix_["X8"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X9"]
        pbm.A[ig,iv] += Float64(2.0)
        iv = ix_["X10"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("C2",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C2")
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(300.0)
        iv = ix_["X3"]
        pbm.A[ig,iv] += Float64(0.09)
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(0.03)
        ig,ig_,_ = s2mpj_ii("C3",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C3")
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(0.326)
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(-101.0)
        iv = ix_["X5"]
        pbm.A[ig,iv] += Float64(200.0)
        iv = ix_["X6"]
        pbm.A[ig,iv] += Float64(0.06)
        iv = ix_["X7"]
        pbm.A[ig,iv] += Float64(0.02)
        ig,ig_,_ = s2mpj_ii("C4",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C4")
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(0.0066667)
        iv = ix_["X3"]
        pbm.A[ig,iv] += Float64(-1.03)
        iv = ix_["X6"]
        pbm.A[ig,iv] += Float64(200.0)
        iv = ix_["X8"]
        pbm.A[ig,iv] += Float64(0.06)
        iv = ix_["X9"]
        pbm.A[ig,iv] += Float64(0.02)
        ig,ig_,_ = s2mpj_ii("C5",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C5")
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(6.6667e-4)
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(-1.01)
        iv = ix_["X7"]
        pbm.A[ig,iv] += Float64(200.0)
        iv = ix_["X9"]
        pbm.A[ig,iv] += Float64(0.06)
        iv = ix_["X10"]
        pbm.A[ig,iv] += Float64(0.02)
        ig,ig_,_ = s2mpj_ii("C6",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C6")
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(0.978)
        iv = ix_["X5"]
        pbm.A[ig,iv] += Float64(-201.0)
        iv = ix_["X11"]
        pbm.A[ig,iv] += Float64(100.0)
        iv = ix_["X12"]
        pbm.A[ig,iv] += Float64(0.03)
        iv = ix_["X13"]
        pbm.A[ig,iv] += Float64(0.01)
        ig,ig_,_ = s2mpj_ii("C7",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C7")
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(0.01)
        iv = ix_["X3"]
        pbm.A[ig,iv] += Float64(0.489)
        iv = ix_["X6"]
        pbm.A[ig,iv] += Float64(-101.03)
        iv = ix_["X12"]
        pbm.A[ig,iv] += Float64(100.0)
        iv = ix_["X14"]
        pbm.A[ig,iv] += Float64(0.03)
        iv = ix_["X15"]
        pbm.A[ig,iv] += Float64(0.01)
        ig,ig_,_ = s2mpj_ii("C8",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C8")
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(0.001)
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(0.489)
        iv = ix_["X7"]
        pbm.A[ig,iv] += Float64(-101.03)
        iv = ix_["X13"]
        pbm.A[ig,iv] += Float64(100.0)
        iv = ix_["X15"]
        pbm.A[ig,iv] += Float64(0.03)
        iv = ix_["X16"]
        pbm.A[ig,iv] += Float64(0.01)
        ig,ig_,_ = s2mpj_ii("C9",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C9")
        iv = ix_["X3"]
        pbm.A[ig,iv] += Float64(0.001)
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(0.01)
        iv = ix_["X9"]
        pbm.A[ig,iv] += Float64(-1.04)
        iv = ix_["X15"]
        pbm.A[ig,iv] += Float64(100.0)
        iv = ix_["X18"]
        pbm.A[ig,iv] += Float64(0.03)
        iv = ix_["X19"]
        pbm.A[ig,iv] += Float64(0.01)
        ig,ig_,_ = s2mpj_ii("C10",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C10")
        iv = ix_["X3"]
        pbm.A[ig,iv] += Float64(0.02)
        iv = ix_["X8"]
        pbm.A[ig,iv] += Float64(-1.06)
        iv = ix_["X14"]
        pbm.A[ig,iv] += Float64(100.0)
        iv = ix_["X17"]
        pbm.A[ig,iv] += Float64(0.03)
        iv = ix_["X19"]
        pbm.A[ig,iv] += Float64(0.01)
        ig,ig_,_ = s2mpj_ii("C11",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C11")
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(0.002)
        iv = ix_["X10"]
        pbm.A[ig,iv] += Float64(-1.02)
        iv = ix_["X16"]
        pbm.A[ig,iv] += Float64(100.0)
        iv = ix_["X19"]
        pbm.A[ig,iv] += Float64(0.03)
        iv = ix_["X20"]
        pbm.A[ig,iv] += Float64(0.01)
        ig,ig_,_ = s2mpj_ii("C12",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C12")
        iv = ix_["X11"]
        pbm.A[ig,iv] += Float64(-2.5742e-6)
        iv = ix_["X13"]
        pbm.A[ig,iv] += Float64(0.00252)
        iv = ix_["X16"]
        pbm.A[ig,iv] += Float64(-0.61975)
        iv = ix_["X20"]
        pbm.A[ig,iv] += Float64(1.03)
        ig,ig_,_ = s2mpj_ii("C13",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C13")
        iv = ix_["X11"]
        pbm.A[ig,iv] += Float64(-0.00257)
        iv = ix_["X12"]
        pbm.A[ig,iv] += Float64(0.25221)
        iv = ix_["X14"]
        pbm.A[ig,iv] += Float64(-6.2)
        iv = ix_["X17"]
        pbm.A[ig,iv] += Float64(1.09)
        ig,ig_,_ = s2mpj_ii("C14",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C14")
        iv = ix_["X11"]
        pbm.A[ig,iv] += Float64(0.00629)
        iv = ix_["X12"]
        pbm.A[ig,iv] += Float64(-0.20555)
        iv = ix_["X13"]
        pbm.A[ig,iv] += Float64(-4.1106)
        iv = ix_["X15"]
        pbm.A[ig,iv] += Float64(101.04)
        iv = ix_["X16"]
        pbm.A[ig,iv] += Float64(505.1)
        iv = ix_["X19"]
        pbm.A[ig,iv] += Float64(-256.72)
        ig,ig_,_ = s2mpj_ii("C15",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C15")
        iv = ix_["X12"]
        pbm.A[ig,iv] += Float64(0.00841)
        iv = ix_["X13"]
        pbm.A[ig,iv] += Float64(-0.08406)
        iv = ix_["X14"]
        pbm.A[ig,iv] += Float64(-0.20667)
        iv = ix_["X16"]
        pbm.A[ig,iv] += Float64(20.658)
        iv = ix_["X18"]
        pbm.A[ig,iv] += Float64(1.07)
        iv = ix_["X19"]
        pbm.A[ig,iv] += Float64(-10.5)
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
        pbm.gconst[ig_["C1"]] = Float64(0.70785)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xupper = fill(1.0,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(1.0),pb.n)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN               3.06435
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
        pb.lincons   = collect(1:length(pbm.congrps))
        pb.pbclass = "C-LLR2-AN-20-15"
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

