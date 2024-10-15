function YAO(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
#    A linear least-sqaure problem with k-convex constraints
#       min (1/2) || f(t) - x ||^2
#    subject to the constraints
#       _ 
#       \/_k  x  >=  0,
#    where  f(t) and  x  are vectors in (n+k)-dimensional space.
# 
#    We choose f(t) = sin(t), x(1) >= 0.08 and fix x(n+i) = 0
# 
#    SIF input: Aixiang Yao, Virginia Tech., May 1995
#               modifications by Nick Gould
# 
#    classification = "C-QLR2-AN-V-V"
# 
#   Number of discretization points
# 
#       Alternative values for the SIF file parameters:
# IE P                   20             $-PARAMETER
# IE P                   200            $-PARAMETER
# IE P                   2000           $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "YAO"

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
            v_["P"] = Int64(20);  #  SIF file default value
        else
            v_["P"] = Int64(args[1]);
        end
# IE k                   2              $-PARAMETER
        if nargin<2
            v_["k"] = Int64(2);  #  SIF file default value
        else
            v_["k"] = Int64(args[2]);
        end
# IE k                   3              $-PARAMETER
# IE k                   4              $-PARAMETER
        v_["1"] = 1
        v_["2"] = 2
        v_["P+1"] = v_["P"]+v_["1"]
        v_["P+k"] = v_["P"]+v_["k"]
        v_["RP"] = Float64(v_["P+k"])
        v_["OVP"] = 1.0/v_["RP"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for i = Int64(v_["1"]):Int64(v_["P+k"])
            iv,ix_,_ = s2mpj_ii("X"*string(i),ix_)
            arrset(pb.xnames,iv,"X"*string(i))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for i = Int64(v_["1"]):Int64(v_["P+k"])
            ig,ig_,_ = s2mpj_ii("S"*string(i),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(i)]
            pbm.A[ig,iv] += Float64(1.0)
            arrset(pbm.gscale,ig,Float64(2.0))
        end
        for i = Int64(v_["1"]):Int64(v_["P"])
            v_["i+1"] = 1+i
            ig,ig_,_ = s2mpj_ii("B"*string(i),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"B"*string(i))
            iv = ix_["X"*string(i)]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["X"*string(Int64(v_["i+1"]))]
            pbm.A[ig,iv] += Float64(-2.0)
            v_["i+2"] = 2+i
            iv = ix_["X"*string(Int64(v_["i+2"]))]
            pbm.A[ig,iv] += Float64(1.0)
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
        for i = Int64(v_["1"]):Int64(v_["P+k"])
            v_["Ri"] = Float64(i)
            v_["iOVP"] = v_["Ri"]*v_["OVP"]
            v_["SINI"] = sin(v_["iOVP"])
            pbm.gconst[ig_["S"*string(i)]] = Float64(v_["SINI"])
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        pb.xlower[ix_["X"*string(Int64(v_["1"]))]] = 0.08
        for i = Int64(v_["P+1"]):Int64(v_["P+k"])
            pb.xlower[ix_["X"*string(i)]] = 0.0
            pb.xupper[ix_["X"*string(i)]] = 0.0
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gSQ",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for i = Int64(v_["1"]):Int64(v_["P+k"])
            ig = ig_["S"*string(i)]
            arrset(pbm.grftype,ig,"gSQ")
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# XL SOLUTION             2.39883D+00   $ (p=20)
# XL SOLUTION             2.01517D+01   $ (p=200)
# XL SOLUTION             1.97705D+02   $ (p=2000)
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
        pb.pbclass = "C-QLR2-AN-V-V"
        pb.x0          = zeros(Float64,pb.n)
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm

# ********************
#  SET UP THE GROUPS *
#  ROUTINE           *
# ********************

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    elseif action == "gSQ"

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

