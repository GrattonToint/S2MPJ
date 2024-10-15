function READING2(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : READING2
#    *********
# 
#    A linear optimal control problem from Nancy Nichols
#    with a given initial condition.
#    This problem arises in tide modelling.
# 
#    Source:
#    S. Lyle and N.K. Nichols,
#    "Numerical Methods for Optimal Control Problems with State Constraints",
#    Numerical Analysis Report 8/91, Dept of Mathematics, 
#    University of Reading, UK.
# 
#    SIF input: Nick Gould, July 1991.
# 
#    classification = "C-LLR2-MN-V-V"
# 
#    Number of discretized points in [0,1]
# 
#       Alternative values for the SIF file parameters:
# IE N                   2              $-PARAMETER
# IE N                   50             $-PARAMETER
# IE N                   100            $-PARAMETER     original value
# IE N                   500            $-PARAMETER
# IE N                   1000           $-PARAMETER
# IE N                   2000           $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "READING2"

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
            v_["N"] = Int64(10);  #  SIF file default value
        else
            v_["N"] = Int64(args[1]);
        end
# IE N                   5000           $-PARAMETER
        v_["PI"] = 3.1415926535
        v_["2PI"] = 2.0*v_["PI"]
        v_["PI**2"] = v_["PI"]*v_["PI"]
        v_["8PI**2"] = 8.0*v_["PI**2"]
        v_["1/8PI**2"] = 1.0/v_["8PI**2"]
        v_["A"] = 0.07716
        v_["1/A"] = 1.0/v_["A"]
        v_["1/2A"] = 2.0*v_["1/A"]
        v_["2A"] = 2.0*v_["A"]
        v_["-2A"] = -1.0*v_["2A"]
        v_["-1/2A"] = 1.0/v_["-2A"]
        v_["N-1"] = -1+v_["N"]
        v_["RN"] = Float64(v_["N"])
        v_["H"] = 1.0/v_["RN"]
        v_["2/H"] = 2.0*v_["RN"]
        v_["H/2"] = 0.5*v_["H"]
        v_["-H/2"] = -1.0*v_["H/2"]
        v_["1/H"] = 1.0*v_["RN"]
        v_["-1/H"] = -1.0*v_["RN"]
        v_["H/8PI**2"] = v_["1/8PI**2"]*v_["H"]
        v_["0"] = 0
        v_["1"] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["0"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("X1u"*string(I),ix_)
            arrset(pb.xnames,iv,"X1u"*string(I))
            iv,ix_,_ = s2mpj_ii("X2u"*string(I),ix_)
            arrset(pb.xnames,iv,"X2u"*string(I))
            iv,ix_,_ = s2mpj_ii("U"*string(I),ix_)
            arrset(pb.xnames,iv,"U"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            v_["RI"] = Float64(I)
            v_["TI"] = v_["RI"]*v_["H"]
            v_["2PITI"] = v_["2PI"]*v_["TI"]
            v_["CTI"] = cos(v_["2PITI"])
            v_["-CCTI"] = v_["CTI"]*v_["-H/2"]
            v_["I-1"] = -1+I
            v_["RI-1"] = Float64(v_["I-1"])
            v_["TI-1"] = v_["RI-1"]*v_["H"]
            v_["2PITI-1"] = v_["2PI"]*v_["TI-1"]
            v_["CTI-1"] = cos(v_["2PITI-1"])
            v_["-CCTI-1"] = v_["CTI-1"]*v_["-H/2"]
            ig,ig_,_ = s2mpj_ii("COST",ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X1u"*string(I)]
            pbm.A[ig,iv] += Float64(v_["-CCTI"])
            iv = ix_["X1u"*string(Int64(v_["I-1"]))]
            pbm.A[ig,iv] += Float64(v_["-CCTI-1"])
            iv = ix_["U"*string(I)]
            pbm.A[ig,iv] += Float64(v_["H/8PI**2"])
            iv = ix_["U"*string(Int64(v_["I-1"]))]
            pbm.A[ig,iv] += Float64(v_["H/8PI**2"])
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            v_["I-1"] = -1+I
            ig,ig_,_ = s2mpj_ii("C1u"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"C1u"*string(I))
            iv = ix_["X1u"*string(I)]
            pbm.A[ig,iv] += Float64(v_["1/H"])
            iv = ix_["X1u"*string(Int64(v_["I-1"]))]
            pbm.A[ig,iv] += Float64(v_["-1/H"])
            iv = ix_["X2u"*string(I)]
            pbm.A[ig,iv] += Float64(-0.5)
            iv = ix_["X2u"*string(Int64(v_["I-1"]))]
            pbm.A[ig,iv] += Float64(-0.5)
            ig,ig_,_ = s2mpj_ii("C2u"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"C2u"*string(I))
            iv = ix_["X2u"*string(I)]
            pbm.A[ig,iv] += Float64(v_["1/H"])
            iv = ix_["X2u"*string(Int64(v_["I-1"]))]
            pbm.A[ig,iv] += Float64(v_["-1/H"])
            iv = ix_["U"*string(I)]
            pbm.A[ig,iv] += Float64(-0.5)
            iv = ix_["U"*string(Int64(v_["I-1"]))]
            pbm.A[ig,iv] += Float64(-0.5)
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
        pb.xlower[ix_["X1u"*string(Int64(v_["0"]))]] = 0.0
        pb.xupper[ix_["X1u"*string(Int64(v_["0"]))]] = 0.0
        pb.xlower[ix_["X2u"*string(Int64(v_["0"]))]] = 0.0
        pb.xupper[ix_["X2u"*string(Int64(v_["0"]))]] = 0.0
        for I = Int64(v_["1"]):Int64(v_["N"])
            pb.xlower[ix_["X1u"*string(I)]] = -Inf
            pb.xupper[ix_["X1u"*string(I)]] = +Inf
            pb.xlower[ix_["X2u"*string(I)]] = -0.125
            pb.xupper[ix_["X2u"*string(I)]] = 0.125
        end
        for I = Int64(v_["0"]):Int64(v_["N"])
            pb.xlower[ix_["U"*string(I)]] = -1.0
            pb.xupper[ix_["U"*string(I)]] = 1.0
        end
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
        pb.pbclass = "C-LLR2-MN-V-V"
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

