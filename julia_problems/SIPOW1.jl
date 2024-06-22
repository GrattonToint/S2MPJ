function SIPOW1(action,args...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : SIPOW1
#    *********
# 
#    This is a discretization of a semi-infinite programming problem, of
#    minimizing the variable x_2 within a circle of radius 1. The circle
#    is replaced by a discrete set of equally-spaced supporting tangents.   
# 
#    Source: problem 1 in
#    M. J. D. Powell,
#    "Log barrier methods for semi-infinite programming calculations"
#    Numerical Analysis Report DAMTP 1992/NA11, U. of Cambridge, UK.
# 
#    SIF input: A. R. Conn and Nick Gould, August 1993
# 
#    classification = "LLR2-AN-2-V"
# 
#    Problem variants: they are identified by the values of M
# 
# IE M                   20 
# IE M                   100 
# IE M                   500 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "SIPOW1"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = eval( Meta.parse( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["M"] = 2000
        v_["1"] = 1
        v_["RM"] = Float64(v_["M"])
        v_["1/RM"] = 1.0/v_["RM"]
        v_["ONE"] = 1.0
        v_["PI/4"] = atan(v_["ONE"])
        v_["2PI"] = 8.0*v_["PI/4"]
        v_["2PI/M"] = v_["2PI"]*v_["1/RM"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        iv,ix_,_ = s2mpj_ii("X1",ix_)
        arrset(pb.xnames,iv,"X1")
        iv,ix_,_ = s2mpj_ii("X2",ix_)
        arrset(pb.xnames,iv,"X2")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(1.0)
        for J = Int64(v_["1"]):Int64(v_["M"])
            v_["RJ"] = Float64(J)
            v_["2PIJ/M"] = v_["2PI/M"]*v_["RJ"]
            v_["COS"] = cos(v_["2PIJ/M"])
            v_["SIN"] = sin(v_["2PIJ/M"])
            ig,ig_,_ = s2mpj_ii("C"*string(J),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"C"*string(J))
            iv = ix_["X1"]
            pbm.A[ig,iv] += Float64(v_["COS"])
            iv = ix_["X2"]
            pbm.A[ig,iv] += Float64(v_["SIN"])
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
        pbm.congrps = findall(x->x!="<>",gtype)
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
        pb.xlower = zeros(Float64,pb.n)
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
        pb.pbclass = "LLR2-AN-2-V"
        return pb, pbm

    #%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    elseif action in  ["fx","fgx","fgHx","cx","cJx","cJHx","cIx","cIJx","cIJHx","cIJxv","fHxv","cJxv","Lxy","Lgxy","LgHxy","LIxy","LIgxy","LIgHxy","LHxyv","LIHxyv"]

        pbm = args[1]
        if pbm.name == name
            pbm.has_globs = [0,0]
            return s2mpj_eval(action,args...)
        else
            println("ERROR: please run "*name*" with action = setup")
            return ntuple(i->undef,args[end])
        end

    else
        println("ERROR: unknown action "*action*" requested from "*name*"%s.jl")
        return ntuple(i->undef,args[end])
    end

end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

