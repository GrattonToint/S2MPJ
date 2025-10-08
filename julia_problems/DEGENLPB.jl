function DEGENLPB(action::String,args::Union{Any}...)
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
#    classification = "C-CLLR2-AN-20-15"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 8 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "DEGENLPB"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling DEGENLPB.")
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
        v_["N"] = 20
        v_["M"] = 15
        v_["1"] = 1
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
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X2"])
        push!(valA,Float64(-0.01))
        push!(irA,ig)
        push!(icA,ix_["X3"])
        push!(valA,Float64(-33.333))
        push!(irA,ig)
        push!(icA,ix_["X4"])
        push!(valA,Float64(-100.0))
        push!(irA,ig)
        push!(icA,ix_["X5"])
        push!(valA,Float64(-0.01))
        push!(irA,ig)
        push!(icA,ix_["X6"])
        push!(valA,Float64(-33.343))
        push!(irA,ig)
        push!(icA,ix_["X7"])
        push!(valA,Float64(-100.01))
        push!(irA,ig)
        push!(icA,ix_["X8"])
        push!(valA,Float64(-33.333))
        push!(irA,ig)
        push!(icA,ix_["X9"])
        push!(valA,Float64(-133.33))
        push!(irA,ig)
        push!(icA,ix_["X10"])
        push!(valA,Float64(-100.0))
        ig,ig_,_ = s2mpj_ii("C1",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C1")
        push!(irA,ig)
        push!(icA,ix_["X1"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X2"])
        push!(valA,Float64(2.0))
        push!(irA,ig)
        push!(icA,ix_["X3"])
        push!(valA,Float64(2.0))
        push!(irA,ig)
        push!(icA,ix_["X4"])
        push!(valA,Float64(2.0))
        push!(irA,ig)
        push!(icA,ix_["X5"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X6"])
        push!(valA,Float64(2.0))
        push!(irA,ig)
        push!(icA,ix_["X7"])
        push!(valA,Float64(2.0))
        push!(irA,ig)
        push!(icA,ix_["X8"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X9"])
        push!(valA,Float64(2.0))
        push!(irA,ig)
        push!(icA,ix_["X10"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("C2",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C2")
        push!(irA,ig)
        push!(icA,ix_["X1"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X2"])
        push!(valA,Float64(300.0))
        push!(irA,ig)
        push!(icA,ix_["X3"])
        push!(valA,Float64(0.09))
        push!(irA,ig)
        push!(icA,ix_["X4"])
        push!(valA,Float64(0.03))
        ig,ig_,_ = s2mpj_ii("C3",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C3")
        push!(irA,ig)
        push!(icA,ix_["X1"])
        push!(valA,Float64(0.326))
        push!(irA,ig)
        push!(icA,ix_["X2"])
        push!(valA,Float64(-101.0))
        push!(irA,ig)
        push!(icA,ix_["X5"])
        push!(valA,Float64(200.0))
        push!(irA,ig)
        push!(icA,ix_["X6"])
        push!(valA,Float64(0.06))
        push!(irA,ig)
        push!(icA,ix_["X7"])
        push!(valA,Float64(0.02))
        ig,ig_,_ = s2mpj_ii("C4",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C4")
        push!(irA,ig)
        push!(icA,ix_["X1"])
        push!(valA,Float64(0.0066667))
        push!(irA,ig)
        push!(icA,ix_["X3"])
        push!(valA,Float64(-1.03))
        push!(irA,ig)
        push!(icA,ix_["X6"])
        push!(valA,Float64(200.0))
        push!(irA,ig)
        push!(icA,ix_["X8"])
        push!(valA,Float64(0.06))
        push!(irA,ig)
        push!(icA,ix_["X9"])
        push!(valA,Float64(0.02))
        ig,ig_,_ = s2mpj_ii("C5",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C5")
        push!(irA,ig)
        push!(icA,ix_["X1"])
        push!(valA,Float64(6.6667e-4))
        push!(irA,ig)
        push!(icA,ix_["X4"])
        push!(valA,Float64(-1.01))
        push!(irA,ig)
        push!(icA,ix_["X7"])
        push!(valA,Float64(200.0))
        push!(irA,ig)
        push!(icA,ix_["X9"])
        push!(valA,Float64(0.06))
        push!(irA,ig)
        push!(icA,ix_["X10"])
        push!(valA,Float64(0.02))
        ig,ig_,_ = s2mpj_ii("C6",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C6")
        push!(irA,ig)
        push!(icA,ix_["X2"])
        push!(valA,Float64(0.978))
        push!(irA,ig)
        push!(icA,ix_["X5"])
        push!(valA,Float64(-201.0))
        push!(irA,ig)
        push!(icA,ix_["X11"])
        push!(valA,Float64(100.0))
        push!(irA,ig)
        push!(icA,ix_["X12"])
        push!(valA,Float64(0.03))
        push!(irA,ig)
        push!(icA,ix_["X13"])
        push!(valA,Float64(0.01))
        ig,ig_,_ = s2mpj_ii("C7",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C7")
        push!(irA,ig)
        push!(icA,ix_["X2"])
        push!(valA,Float64(0.01))
        push!(irA,ig)
        push!(icA,ix_["X3"])
        push!(valA,Float64(0.489))
        push!(irA,ig)
        push!(icA,ix_["X6"])
        push!(valA,Float64(-101.03))
        push!(irA,ig)
        push!(icA,ix_["X12"])
        push!(valA,Float64(100.0))
        push!(irA,ig)
        push!(icA,ix_["X14"])
        push!(valA,Float64(0.03))
        push!(irA,ig)
        push!(icA,ix_["X15"])
        push!(valA,Float64(0.01))
        ig,ig_,_ = s2mpj_ii("C8",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C8")
        push!(irA,ig)
        push!(icA,ix_["X2"])
        push!(valA,Float64(0.001))
        push!(irA,ig)
        push!(icA,ix_["X4"])
        push!(valA,Float64(0.489))
        push!(irA,ig)
        push!(icA,ix_["X7"])
        push!(valA,Float64(-101.03))
        push!(irA,ig)
        push!(icA,ix_["X13"])
        push!(valA,Float64(100.0))
        push!(irA,ig)
        push!(icA,ix_["X15"])
        push!(valA,Float64(0.03))
        push!(irA,ig)
        push!(icA,ix_["X16"])
        push!(valA,Float64(0.01))
        ig,ig_,_ = s2mpj_ii("C9",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C9")
        push!(irA,ig)
        push!(icA,ix_["X3"])
        push!(valA,Float64(0.001))
        push!(irA,ig)
        push!(icA,ix_["X4"])
        push!(valA,Float64(0.01))
        push!(irA,ig)
        push!(icA,ix_["X9"])
        push!(valA,Float64(-1.04))
        push!(irA,ig)
        push!(icA,ix_["X15"])
        push!(valA,Float64(100.0))
        push!(irA,ig)
        push!(icA,ix_["X18"])
        push!(valA,Float64(0.03))
        push!(irA,ig)
        push!(icA,ix_["X19"])
        push!(valA,Float64(0.01))
        ig,ig_,_ = s2mpj_ii("C10",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C10")
        push!(irA,ig)
        push!(icA,ix_["X3"])
        push!(valA,Float64(0.02))
        push!(irA,ig)
        push!(icA,ix_["X8"])
        push!(valA,Float64(-1.06))
        push!(irA,ig)
        push!(icA,ix_["X14"])
        push!(valA,Float64(100.0))
        push!(irA,ig)
        push!(icA,ix_["X17"])
        push!(valA,Float64(0.03))
        push!(irA,ig)
        push!(icA,ix_["X19"])
        push!(valA,Float64(0.01))
        ig,ig_,_ = s2mpj_ii("C11",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C11")
        push!(irA,ig)
        push!(icA,ix_["X4"])
        push!(valA,Float64(0.002))
        push!(irA,ig)
        push!(icA,ix_["X10"])
        push!(valA,Float64(-1.02))
        push!(irA,ig)
        push!(icA,ix_["X16"])
        push!(valA,Float64(100.0))
        push!(irA,ig)
        push!(icA,ix_["X19"])
        push!(valA,Float64(0.03))
        push!(irA,ig)
        push!(icA,ix_["X20"])
        push!(valA,Float64(0.01))
        ig,ig_,_ = s2mpj_ii("C12",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C12")
        push!(irA,ig)
        push!(icA,ix_["X11"])
        push!(valA,Float64(-2.5742e-6))
        push!(irA,ig)
        push!(icA,ix_["X13"])
        push!(valA,Float64(0.00252))
        push!(irA,ig)
        push!(icA,ix_["X16"])
        push!(valA,Float64(-0.61975))
        push!(irA,ig)
        push!(icA,ix_["X20"])
        push!(valA,Float64(1.03))
        ig,ig_,_ = s2mpj_ii("C13",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C13")
        push!(irA,ig)
        push!(icA,ix_["X11"])
        push!(valA,Float64(-0.00257))
        push!(irA,ig)
        push!(icA,ix_["X12"])
        push!(valA,Float64(0.25221))
        push!(irA,ig)
        push!(icA,ix_["X14"])
        push!(valA,Float64(-6.2))
        push!(irA,ig)
        push!(icA,ix_["X17"])
        push!(valA,Float64(1.09))
        ig,ig_,_ = s2mpj_ii("C14",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C14")
        push!(irA,ig)
        push!(icA,ix_["X11"])
        push!(valA,Float64(0.00629))
        push!(irA,ig)
        push!(icA,ix_["X12"])
        push!(valA,Float64(-0.20555))
        push!(irA,ig)
        push!(icA,ix_["X13"])
        push!(valA,Float64(-4.1106))
        push!(irA,ig)
        push!(icA,ix_["X15"])
        push!(valA,Float64(101.04))
        push!(irA,ig)
        push!(icA,ix_["X16"])
        push!(valA,Float64(505.1))
        push!(irA,ig)
        push!(icA,ix_["X19"])
        push!(valA,Float64(-256.72))
        ig,ig_,_ = s2mpj_ii("C15",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C15")
        push!(irA,ig)
        push!(icA,ix_["X12"])
        push!(valA,Float64(0.00841))
        push!(irA,ig)
        push!(icA,ix_["X13"])
        push!(valA,Float64(-0.08406))
        push!(irA,ig)
        push!(icA,ix_["X14"])
        push!(valA,Float64(-0.20667))
        push!(irA,ig)
        push!(icA,ix_["X16"])
        push!(valA,Float64(20.658))
        push!(irA,ig)
        push!(icA,ix_["X18"])
        push!(valA,Float64(1.07))
        push!(irA,ig)
        push!(icA,ix_["X19"])
        push!(valA,Float64(-10.5))
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
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = collect(1:length(pbm.congrps))
        pb.pbclass = "C-CLLR2-AN-20-15"
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

