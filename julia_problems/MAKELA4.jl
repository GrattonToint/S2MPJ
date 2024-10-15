function MAKELA4(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : MAKELA4
#    *********
# 
#    A nonlinear minmax problem in twenty variables.
# 
#    Source: 
#    M.M. Makela,
#    "Nonsmooth optimization",
#    Ph.D. thesis, Jyvaskyla University, 1990
# 
#    SIF input: Ph. Toint, Nov 1993.
# 
#    classification = "C-LLR2-AN-21-40"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "MAKELA4"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["1"] = 1
        v_["20"] = 20
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["20"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        iv,ix_,_ = s2mpj_ii("U",ix_)
        arrset(pb.xnames,iv,"U")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["U"]
        pbm.A[ig,iv] += Float64(1.0)
        for I = Int64(v_["1"]):Int64(v_["20"])
            ig,ig_,_ = s2mpj_ii("F"*string(I),ig_)
            arrset(gtype,ig,"<=")
            arrset(pb.cnames,ig,"F"*string(I))
            iv = ix_["U"]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("MF"*string(I),ig_)
            arrset(gtype,ig,"<=")
            arrset(pb.cnames,ig,"MF"*string(I))
            iv = ix_["U"]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
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
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        if haskey(ix_,"X1")
            pb.x0[ix_["X1"]] = Float64(1.0)
        else
            pb.y0[findfirst(x->x==ig_["X1"],pbm.congrps)] = Float64(1.0)
        end
        if haskey(ix_,"X2")
            pb.x0[ix_["X2"]] = Float64(2.0)
        else
            pb.y0[findfirst(x->x==ig_["X2"],pbm.congrps)] = Float64(2.0)
        end
        if haskey(ix_,"X3")
            pb.x0[ix_["X3"]] = Float64(3.0)
        else
            pb.y0[findfirst(x->x==ig_["X3"],pbm.congrps)] = Float64(3.0)
        end
        if haskey(ix_,"X4")
            pb.x0[ix_["X4"]] = Float64(4.0)
        else
            pb.y0[findfirst(x->x==ig_["X4"],pbm.congrps)] = Float64(4.0)
        end
        if haskey(ix_,"X5")
            pb.x0[ix_["X5"]] = Float64(5.0)
        else
            pb.y0[findfirst(x->x==ig_["X5"],pbm.congrps)] = Float64(5.0)
        end
        if haskey(ix_,"X6")
            pb.x0[ix_["X6"]] = Float64(6.0)
        else
            pb.y0[findfirst(x->x==ig_["X6"],pbm.congrps)] = Float64(6.0)
        end
        if haskey(ix_,"X7")
            pb.x0[ix_["X7"]] = Float64(7.0)
        else
            pb.y0[findfirst(x->x==ig_["X7"],pbm.congrps)] = Float64(7.0)
        end
        if haskey(ix_,"X8")
            pb.x0[ix_["X8"]] = Float64(8.0)
        else
            pb.y0[findfirst(x->x==ig_["X8"],pbm.congrps)] = Float64(8.0)
        end
        if haskey(ix_,"X9")
            pb.x0[ix_["X9"]] = Float64(9.0)
        else
            pb.y0[findfirst(x->x==ig_["X9"],pbm.congrps)] = Float64(9.0)
        end
        if haskey(ix_,"X10")
            pb.x0[ix_["X10"]] = Float64(10.0)
        else
            pb.y0[findfirst(x->x==ig_["X10"],pbm.congrps)] = Float64(10.0)
        end
        if haskey(ix_,"X11")
            pb.x0[ix_["X11"]] = Float64(-11.0)
        else
            pb.y0[findfirst(x->x==ig_["X11"],pbm.congrps)] = Float64(-11.0)
        end
        if haskey(ix_,"X12")
            pb.x0[ix_["X12"]] = Float64(-12.0)
        else
            pb.y0[findfirst(x->x==ig_["X12"],pbm.congrps)] = Float64(-12.0)
        end
        if haskey(ix_,"X13")
            pb.x0[ix_["X13"]] = Float64(-13.0)
        else
            pb.y0[findfirst(x->x==ig_["X13"],pbm.congrps)] = Float64(-13.0)
        end
        if haskey(ix_,"X14")
            pb.x0[ix_["X14"]] = Float64(-14.0)
        else
            pb.y0[findfirst(x->x==ig_["X14"],pbm.congrps)] = Float64(-14.0)
        end
        if haskey(ix_,"X15")
            pb.x0[ix_["X15"]] = Float64(-15.0)
        else
            pb.y0[findfirst(x->x==ig_["X15"],pbm.congrps)] = Float64(-15.0)
        end
        if haskey(ix_,"X16")
            pb.x0[ix_["X16"]] = Float64(-16.0)
        else
            pb.y0[findfirst(x->x==ig_["X16"],pbm.congrps)] = Float64(-16.0)
        end
        if haskey(ix_,"X17")
            pb.x0[ix_["X17"]] = Float64(-17.0)
        else
            pb.y0[findfirst(x->x==ig_["X17"],pbm.congrps)] = Float64(-17.0)
        end
        if haskey(ix_,"X18")
            pb.x0[ix_["X18"]] = Float64(-18.0)
        else
            pb.y0[findfirst(x->x==ig_["X18"],pbm.congrps)] = Float64(-18.0)
        end
        if haskey(ix_,"X19")
            pb.x0[ix_["X19"]] = Float64(-19.0)
        else
            pb.y0[findfirst(x->x==ig_["X19"],pbm.congrps)] = Float64(-19.0)
        end
        if haskey(ix_,"X20")
            pb.x0[ix_["X20"]] = Float64(-20.0)
        else
            pb.y0[findfirst(x->x==ig_["X20"],pbm.congrps)] = Float64(-20.0)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = collect(1:length(pbm.congrps))
        pb.pbclass = "C-LLR2-AN-21-40"
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

