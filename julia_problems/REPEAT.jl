function REPEAT(action,args...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : REPEAT
#    *********
# 
#    This problem is to find the nearest feasible point to 2n+1 inconsistent
#    linear equations subject to bounds
# 
#    Source: blue-cheese delerium
# 
#    SIF input: Nick Gould, December 2020.
# 
#    classification = "NLR2-AN-V-V"
# 
#    N is the number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "REPEAT"

    if action == "setup"
        pbm          = PBM(name)
        pb           = PB(name)
        pb.sifpbname = "REPEAT"
        nargin       = length(args)
        pbm.call     = eval( Meta.parse( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        if nargin<1
            v_["N"] = 10;  #  SIF file default value
        else
            v_["N"] = args[1];
        end
        v_["N-1"] = -1+v_["N"]
        v_["1"] = 1
        v_["2"] = 2
        v_["100"] = 100
        v_["N/2"] = trunc(Int,(v_["N"]/v_["2"]))
        v_["N/100"] = trunc(Int,(v_["N"]/v_["100"]))
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        xscale  = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2x_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            v_["I+1"] = 1+I
            ig,ig_,_ = s2x_ii("C"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"C"*string(I))
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += 1.0
            iv = ix_["X"*string(Int64(v_["I+1"]))]
            pbm.A[ig,iv] += 1.0
        end
        ig,ig_,_ = s2x_ii("C"*string(Int64(v_["N"])),ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C"*string(Int64(v_["N"])))
        iv = ix_["X"*string(Int64(v_["N"]))]
        pbm.A[ig,iv] += 1.0
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            v_["I+1"] = 1+I
            ig,ig_,_ = s2x_ii("R"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"R"*string(I))
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += 1.0
            iv = ix_["X"*string(Int64(v_["I+1"]))]
            pbm.A[ig,iv] += 1.0
        end
        ig,ig_,_ = s2x_ii("R"*string(Int64(v_["N"])),ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"R"*string(Int64(v_["N"])))
        iv = ix_["X"*string(Int64(v_["N"]))]
        pbm.A[ig,iv] += 1.0
        for I = Int64(v_["1"]):Int64(v_["N/100"]):Int64(v_["N"])
            v_["RI"] =Float64(I)
            ig,ig_,_ = s2x_ii("E",ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"E")
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += v_["RI"]
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
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            pbm.gconst[ig_["C"*string(I)]] = 2.0
        end
        pbm.gconst[ig_["C"*string(Int64(v_["N"]))]] = 1.0
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            pbm.gconst[ig_["R"*string(I)]] = 4.0
        end
        pbm.gconst[ig_["R"*string(Int64(v_["N"]))]] = 3.0
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = fill(-1.0,pb.n)
        pb.xupper[ix_["X"*string(Int64(v_["1"]))]] = 0.0
        pb.xlower[ix_["X"*string(Int64(v_["2"]))]] = 3.0
        pb.xupper[ix_["X"*string(Int64(v_["N/2"]))]] = 0.0
        pb.xlower[ix_["X"*string(Int64(v_["N-1"]))]] = 3.0
        pb.xupper[ix_["X"*string(Int64(v_["N"]))]] = 0.0
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(0.0,pb.n,)
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
        pb.pbclass = "NLR2-AN-V-V"
        return pb, pbm

    #%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    elseif action in  ["fx","fgx","fgHx","cx","cJx","cJHx","cIx","cIJx","cIJHx","cIJxv","fHxv","cJxv","Lxy","Lgxy","LgHxy","LIxy","LIgxy","LIgHxy","LHxyv","LIHxyv"]

        pbm = args[1]
        if pbm.name == name
            pbm.has_globs = [0,0]
            return s2x_eval(action,args...)
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

