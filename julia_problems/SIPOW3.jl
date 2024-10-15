function SIPOW3(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : SIPOW3
#    *********
# 
#    This is a discretization of a one sided approximation problem of
#    approximating the function xi * xi * eta by a linear polynomial
#    on the boundary of the unit square [0,1]x[0,1].
# 
#    Source: problem 3 in
#    M. J. D. Powell,
#    "Log barrier methods for semi-infinite programming calculations"
#    Numerical Analysis Report DAMTP 1992/NA11, U. of Cambridge, UK.
# 
#    SIF input: A. R. Conn and Nick Gould, August 1993
# 
#    classification = "C-LLR2-AN-4-V"
# 
#    Problem variants: they are identified by the values of M (even)
# 
# IE M                   20 
# IE M                   100 
# IE M                   500 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "SIPOW3"

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
        v_["RM"] = Float64(v_["M"])
        v_["RM/8"] = 0.125*v_["RM"]
        v_["RM/4"] = 0.25*v_["RM"]
        v_["3RM/8"] = 0.375*v_["RM"]
        v_["M/8"] = trunc(Int,v_["RM/8"])
        v_["M/4"] = trunc(Int,v_["RM/4"])
        v_["3M/8"] = trunc(Int,v_["3RM/8"])
        v_["M/8+1"] = 1+v_["M/8"]
        v_["M/4+1"] = 1+v_["M/4"]
        v_["3M/8+1"] = 1+v_["3M/8"]
        v_["M/2"] = trunc(Int,(v_["M"]/v_["2"]))
        v_["M/2+1"] = 1+v_["M/2"]
        v_["RM"] = Float64(v_["M"])
        v_["STEP"] = 8.0/v_["RM"]
        for J = Int64(v_["1"]):Int64(v_["M/2"])
            v_["I"] = -1+J
            v_["RI"] = Float64(v_["I"])
            v_["XI"*string(J)] = v_["RI"]*v_["STEP"]
        end
        for J = Int64(v_["1"]):Int64(v_["M/8"])
            v_["RJ"] = Float64(J)
            v_["ETA"*string(J)] = v_["XI"*string(J)]
            v_["XI"*string(J)] = 0.0
        end
        for J = Int64(v_["M/8+1"]):Int64(v_["M/4"])
            v_["RJ"] = Float64(J)
            v_["XI"*string(J)] = -1.0+v_["XI"*string(J)]
            v_["ETA"*string(J)] = 1.0
        end
        for J = Int64(v_["M/4+1"]):Int64(v_["3M/8"])
            v_["RJ"] = Float64(J)
            v_["ETA"*string(J)] = -2.0+v_["XI"*string(J)]
            v_["XI"*string(J)] = 1.0
        end
        for J = Int64(v_["3M/8+1"]):Int64(v_["M/2"])
            v_["RJ"] = Float64(J)
            v_["XI"*string(J)] = -3.0+v_["XI"*string(J)]
            v_["ETA"*string(J)] = 0.0
        end
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        iv,ix_,_ = s2mpj_ii("X1",ix_)
        arrset(pb.xnames,iv,"X1")
        iv,ix_,_ = s2mpj_ii("X2",ix_)
        arrset(pb.xnames,iv,"X2")
        iv,ix_,_ = s2mpj_ii("X3",ix_)
        arrset(pb.xnames,iv,"X3")
        iv,ix_,_ = s2mpj_ii("X4",ix_)
        arrset(pb.xnames,iv,"X4")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(1.0)
        for J = Int64(v_["1"]):Int64(v_["M/2"])
            ig,ig_,_ = s2mpj_ii("C"*string(J),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"C"*string(J))
            iv = ix_["X1"]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["X4"]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["X2"]
            pbm.A[ig,iv] += Float64(v_["XI"*string(J)])
            iv = ix_["X3"]
            pbm.A[ig,iv] += Float64(v_["ETA"*string(J)])
        end
        for J = Int64(v_["1"]):Int64(v_["M/2"])
            v_["J+"] = v_["M/2"]+J
            ig,ig_,_ = s2mpj_ii("C"*string(Int64(v_["J+"])),ig_)
            arrset(gtype,ig,"<=")
            arrset(pb.cnames,ig,"C"*string(Int64(v_["J+"])))
            iv = ix_["X1"]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("C"*string(Int64(v_["J+"])),ig_)
            arrset(gtype,ig,"<=")
            arrset(pb.cnames,ig,"C"*string(Int64(v_["J+"])))
            iv = ix_["X2"]
            pbm.A[ig,iv] += Float64(v_["XI"*string(J)])
            ig,ig_,_ = s2mpj_ii("C"*string(Int64(v_["J+"])),ig_)
            arrset(gtype,ig,"<=")
            arrset(pb.cnames,ig,"C"*string(Int64(v_["J+"])))
            iv = ix_["X3"]
            pbm.A[ig,iv] += Float64(v_["ETA"*string(J)])
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
        for J = Int64(v_["1"]):Int64(v_["M/2"])
            v_["J+"] = v_["M/2"]+J
            v_["XIXI"] = v_["XI"*string(J)]*v_["XI"*string(J)]
            v_["XIXIETA"] = v_["XIXI"]*v_["ETA"*string(J)]
            pbm.gconst[ig_["C"*string(J)]] = Float64(v_["XIXIETA"])
            pbm.gconst[ig_["C"*string(Int64(v_["J+"]))]] = Float64(v_["XIXIETA"])
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        if haskey(ix_,"X1")
            pb.x0[ix_["X1"]] = Float64(-0.1)
        else
            pb.y0[findfirst(x->x==ig_["X1"],pbm.congrps)] = Float64(-0.1)
        end
        if haskey(ix_,"X2")
            pb.x0[ix_["X2"]] = Float64(0.0)
        else
            pb.y0[findfirst(x->x==ig_["X2"],pbm.congrps)] = Float64(0.0)
        end
        if haskey(ix_,"X3")
            pb.x0[ix_["X3"]] = Float64(0.0)
        else
            pb.y0[findfirst(x->x==ig_["X3"],pbm.congrps)] = Float64(0.0)
        end
        if haskey(ix_,"X4")
            pb.x0[ix_["X4"]] = Float64(1.2)
        else
            pb.y0[findfirst(x->x==ig_["X4"],pbm.congrps)] = Float64(1.2)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
# LO SOLUTION            3.0315716D-1 ! m = 20
# LO SOLUTION            5.0397238D-1 ! m = 100
# LO SOLUTION            5.3016386D-1 ! m = 500
# LO SOLUTION            5.3465470D-1 ! m = 2000
# LO SOLUTION            5.3564207D-1 ! m = 10000
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        pb.clower[pb.nle+pb.neq+1:pb.m] = zeros(Float64,pb.nge)
        pb.cupper[1:pb.nge] = fill(Inf,pb.nge)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = collect(1:length(pbm.congrps))
        pb.pbclass = "C-LLR2-AN-4-V"
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

