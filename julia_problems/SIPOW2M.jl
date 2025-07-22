function SIPOW2M(action::String,args::Union{Any}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : SIPOW2M
#    *********
# 
#    This is a discretization of a semi-infinite programming problem, of
#    minimizing the variable x_2 within a circle of radius 1. The circle
#    is replaced by a discrete set of equally-spaced supporting tangents.   
#    The symmetry in SIPOW1.SIF is imposed by replacing those constraints
#    by an alternative set.
# 
#    A modification proposed by Powell, section 6.
# 
#    Source: problem 2 - modified - in
#    M. J. D. Powell,
#    "Log barrier methods for semi-infinite programming calculations"
#    Numerical Analysis Report DAMTP 1992/NA11, U. of Cambridge, UK.
# 
#    SIF input: A. R. Conn and Nick Gould, August 1993
# 
#    classification = "C-CLLR2-AN-2-V"
# 
#    Problem variants: they are identified by the values of M (even)
# 
# IE M                   20 
# IE M                   100 
# IE M                   500 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 22 VII 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "SIPOW2M"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling SIPOW2M.")
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
        v_["M"] = 2000
        v_["1"] = 1
        v_["2"] = 2
        v_["M/2"] = trunc(Int,(v_["M"]/v_["2"]))
        v_["M/2+1"] = 1+v_["M/2"]
        v_["RM"] = Float64(v_["M"])
        v_["1/RM"] = 1.0/v_["RM"]
        v_["ONE"] = 1.0
        v_["PI/4"] = atan(v_["ONE"])
        v_["4PI"] = 16.0*v_["PI/4"]
        v_["4PI/M"] = v_["4PI"]*v_["1/RM"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        irA   = Int64[]
        icA   = Int64[]
        valA  = Float64[]
        iv,ix_,_ = s2mpj_ii("X1",ix_)
        arrset(pb.xnames,iv,"X1")
        iv,ix_,_ = s2mpj_ii("X2",ix_)
        arrset(pb.xnames,iv,"X2")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X2"])
        push!(valA,Float64(1.0))
        for J = Int64(v_["1"]):Int64(v_["M/2"])
            v_["RJ"] = Float64(J)
            v_["RJ+1/2"] = 0.5+v_["RJ"]
            v_["4PIJ+/M"] = v_["4PI/M"]*v_["RJ+1/2"]
            v_["COS"] = cos(v_["4PIJ+/M"])
            v_["SIN"] = sin(v_["4PIJ+/M"])
            ig,ig_,_ = s2mpj_ii("C"*string(J),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"C"*string(J))
            push!(irA,ig)
            push!(icA,ix_["X1"])
            push!(valA,Float64(v_["COS"]))
            push!(irA,ig)
            push!(icA,ix_["X2"])
            push!(valA,Float64(v_["SIN"]))
        end
        for J = Int64(v_["M/2+1"]):Int64(v_["M"])
            ig,ig_,_ = s2mpj_ii("C"*string(J),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"C"*string(J))
            push!(irA,ig)
            push!(icA,ix_["X1"])
            push!(valA,Float64(1.0))
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
        for J = Int64(v_["1"]):Int64(v_["M"])
            pbm.gconst[ig_["C"*string(J)]] = Float64(-1.0)
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        if haskey(ix_,"X1")
            pb.x0[ix_["X1"]] = Float64(0.8)
        else
            pb.y0[findfirst(x->x==ig_["X1"],pbm.congrps)] = Float64(0.8)
        end
        if haskey(ix_,"X2")
            pb.x0[ix_["X2"]] = Float64(0.5)
        else
            pb.y0[findfirst(x->x==ig_["X2"],pbm.congrps)] = Float64(0.5)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
# LO SOLUTION            -1.0
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+pb.neq+1:pb.m] = zeros(Float64,pb.nge)
        pb.cupper[1:pb.nge] = fill(Inf,pb.nge)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = collect(1:length(pbm.congrps))
        pb.pbclass = "C-CLLR2-AN-2-V"
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

