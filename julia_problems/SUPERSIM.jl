function SUPERSIM(action::String,args::Union{Any}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : SUPERSIM
#    *********
#    A super simple linear program (version 2).
# 
#    Source:
#    A.R. Conn, N. Gould and Ph.L. Toint,
#    "The LANCELOT User's Manual",
#    Dept of Maths, FUNDP, 1991.
# 
#    SIF input: Ph. Toint, Jan 1991.
# 
#    classification = "C-CLLR2-AN-2-2"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 21 VI 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "SUPERSIM"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling SUPERSIM.")
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
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        irA   = Int64[]
        icA   = Int64[]
        valA  = Float64[]
        iv,ix_,_ = s2mpj_ii("x",ix_)
        arrset(pb.xnames,iv,"x")
        iv,ix_,_ = s2mpj_ii("y",ix_)
        arrset(pb.xnames,iv,"y")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype = String[]
        ig,ig_,_ = s2mpj_ii("Object",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["x"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("Cautious",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"Cautious")
        push!(irA,ig)
        push!(icA,ix_["x"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["y"])
        push!(valA,Float64(2.0))
        ig,ig_,_ = s2mpj_ii("Daring",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"Daring")
        push!(irA,ig)
        push!(icA,ix_["x"])
        push!(valA,Float64(2.0))
        push!(irA,ig)
        push!(icA,ix_["y"])
        push!(valA,Float64(1.0))
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
        v_["Level"] = 2.0
        pbm.gconst[ig_["Cautious"]] = Float64(v_["Level"])
        pbm.gconst[ig_["Daring"]] = Float64(v_["Level"])
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower[ix_["y"]] = -Inf
        pb.xupper[ix_["y"]] = +Inf
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
        pb.pbclass = "C-CLLR2-AN-2-2"
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

