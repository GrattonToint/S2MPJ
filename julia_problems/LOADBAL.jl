function LOADBAL(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    The problem arises in the field of computer networks and parallel
#    computation.  It deals with the static load balancing in a tree
#    computer network with two-way traffic.  A set of heterogeneous host
#    computers are interconnected, in which each node processes jobs (the 
#    jobs arriving at each node according to a time invariant Poisson process) 
#    locally or sends it to a remote node,.  In the latter case, there is a
#    communication delay of forwarding the job and getting a response back.
#    The problem is then to minimize the mean response time of a job.
# 
#    The example considered here features 11 computers arranged as follows:
# 
#          1      6      9
#           \     |     /
#            \    |    /
#         2---4---5---8---10
#            /    |    \
#           /     |     \
#          3      7      11
# 
#    Source:
#    J. Li and H. Kameda,
#    "Optimal load balancing in tree network with two-way traffic",
#    Computer networks and ISDN systems, vol. 25, pp. 1335-1348, 1993.
# 
#    SIF input: Masha Sosonkina, Virginia Tech., 1995.
# 
#    classification = "C-OLR2-MN-31-31"
# 
#  Parameter assignment.
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "LOADBAL"

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
        v_["P1"] = 3
        v_["N"] = 11
        v_["NLINK"] = 20
        v_["NLINK-3"] = 17
        v_["NLINK-4"] = 16
        v_["4C"] = 4
        v_["5C"] = 5
        v_["6C"] = 6
        v_["7C"] = 7
        v_["8C"] = 8
        v_["FI"] = 514.0
        v_["0.2*FI"] = 0.2*v_["FI"]
        v_["0.0125*FI"] = 0.0125*v_["FI"]
        v_["0.05*FI"] = 0.05*v_["FI"]
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            ig,ig_,_ = s2mpj_ii("F"*string(I),ig_)
            arrset(gtype,ig,"<>")
            arrset(pbm.gscale,ig,Float64(v_["0.2*FI"]))
        end
        for I = Int64(v_["1"]):Int64(v_["NLINK"])
            ig,ig_,_ = s2mpj_ii("CNST"*string(I),ig_)
            arrset(gtype,ig,"<=")
            arrset(pb.cnames,ig,"CNST"*string(I))
            ig,ig_,_ = s2mpj_ii("GA"*string(I),ig_)
            arrset(gtype,ig,"<>")
            arrset(pbm.gscale,ig,Float64(v_["0.0125*FI"]))
            ig,ig_,_ = s2mpj_ii("GB"*string(I),ig_)
            arrset(gtype,ig,"<>")
            arrset(pbm.gscale,ig,Float64(v_["0.05*FI"]))
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            ig,ig_,_ = s2mpj_ii("N"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"N"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        ngrp   = length(ig_)
        iv,ix_,_ = s2mpj_ii("X4,1",ix_)
        arrset(pb.xnames,iv,"X4,1")
        ig = ig_["N1"]
        pbm.A[ig,iv] += Float64(1.00)
        ig = ig_["N4"]
        pbm.A[ig,iv] += Float64(-1.00)
        iv,ix_,_ = s2mpj_ii("X1,4",ix_)
        arrset(pb.xnames,iv,"X1,4")
        ig = ig_["N1"]
        pbm.A[ig,iv] += Float64(-1.00)
        ig = ig_["N4"]
        pbm.A[ig,iv] += Float64(1.00)
        iv,ix_,_ = s2mpj_ii("X4,2",ix_)
        arrset(pb.xnames,iv,"X4,2")
        ig = ig_["N2"]
        pbm.A[ig,iv] += Float64(1.00)
        ig = ig_["N4"]
        pbm.A[ig,iv] += Float64(-1.00)
        iv,ix_,_ = s2mpj_ii("X2,4",ix_)
        arrset(pb.xnames,iv,"X2,4")
        ig = ig_["N2"]
        pbm.A[ig,iv] += Float64(-1.00)
        ig = ig_["N4"]
        pbm.A[ig,iv] += Float64(1.00)
        iv,ix_,_ = s2mpj_ii("X4,3",ix_)
        arrset(pb.xnames,iv,"X4,3")
        ig = ig_["N3"]
        pbm.A[ig,iv] += Float64(1.00)
        ig = ig_["N4"]
        pbm.A[ig,iv] += Float64(-1.00)
        iv,ix_,_ = s2mpj_ii("X3,4",ix_)
        arrset(pb.xnames,iv,"X3,4")
        ig = ig_["N3"]
        pbm.A[ig,iv] += Float64(-1.00)
        ig = ig_["N4"]
        pbm.A[ig,iv] += Float64(1.00)
        iv,ix_,_ = s2mpj_ii("X4,5",ix_)
        arrset(pb.xnames,iv,"X4,5")
        ig = ig_["N5"]
        pbm.A[ig,iv] += Float64(1.00)
        ig = ig_["N4"]
        pbm.A[ig,iv] += Float64(-1.00)
        iv,ix_,_ = s2mpj_ii("X5,4",ix_)
        arrset(pb.xnames,iv,"X5,4")
        ig = ig_["N5"]
        pbm.A[ig,iv] += Float64(-1.00)
        ig = ig_["N4"]
        pbm.A[ig,iv] += Float64(1.00)
        iv,ix_,_ = s2mpj_ii("X5,6",ix_)
        arrset(pb.xnames,iv,"X5,6")
        ig = ig_["N6"]
        pbm.A[ig,iv] += Float64(1.00)
        ig = ig_["N5"]
        pbm.A[ig,iv] += Float64(-1.00)
        iv,ix_,_ = s2mpj_ii("X6,5",ix_)
        arrset(pb.xnames,iv,"X6,5")
        ig = ig_["N6"]
        pbm.A[ig,iv] += Float64(-1.00)
        ig = ig_["N5"]
        pbm.A[ig,iv] += Float64(1.00)
        iv,ix_,_ = s2mpj_ii("X5,7",ix_)
        arrset(pb.xnames,iv,"X5,7")
        ig = ig_["N7"]
        pbm.A[ig,iv] += Float64(1.00)
        ig = ig_["N5"]
        pbm.A[ig,iv] += Float64(-1.00)
        iv,ix_,_ = s2mpj_ii("X7,5",ix_)
        arrset(pb.xnames,iv,"X7,5")
        ig = ig_["N7"]
        pbm.A[ig,iv] += Float64(-1.00)
        ig = ig_["N5"]
        pbm.A[ig,iv] += Float64(1.00)
        iv,ix_,_ = s2mpj_ii("X5,8",ix_)
        arrset(pb.xnames,iv,"X5,8")
        ig = ig_["N8"]
        pbm.A[ig,iv] += Float64(1.00)
        ig = ig_["N5"]
        pbm.A[ig,iv] += Float64(-1.00)
        iv,ix_,_ = s2mpj_ii("X8,5",ix_)
        arrset(pb.xnames,iv,"X8,5")
        ig = ig_["N8"]
        pbm.A[ig,iv] += Float64(-1.00)
        ig = ig_["N5"]
        pbm.A[ig,iv] += Float64(1.00)
        iv,ix_,_ = s2mpj_ii("X8,9",ix_)
        arrset(pb.xnames,iv,"X8,9")
        ig = ig_["N9"]
        pbm.A[ig,iv] += Float64(1.00)
        ig = ig_["N8"]
        pbm.A[ig,iv] += Float64(-1.00)
        iv,ix_,_ = s2mpj_ii("X9,8",ix_)
        arrset(pb.xnames,iv,"X9,8")
        ig = ig_["N9"]
        pbm.A[ig,iv] += Float64(-1.00)
        ig = ig_["N8"]
        pbm.A[ig,iv] += Float64(1.00)
        iv,ix_,_ = s2mpj_ii("X8,10",ix_)
        arrset(pb.xnames,iv,"X8,10")
        ig = ig_["N10"]
        pbm.A[ig,iv] += Float64(1.00)
        ig = ig_["N8"]
        pbm.A[ig,iv] += Float64(-1.00)
        iv,ix_,_ = s2mpj_ii("X10,8",ix_)
        arrset(pb.xnames,iv,"X10,8")
        ig = ig_["N10"]
        pbm.A[ig,iv] += Float64(-1.00)
        ig = ig_["N8"]
        pbm.A[ig,iv] += Float64(1.00)
        iv,ix_,_ = s2mpj_ii("X8,11",ix_)
        arrset(pb.xnames,iv,"X8,11")
        ig = ig_["N11"]
        pbm.A[ig,iv] += Float64(1.00)
        ig = ig_["N8"]
        pbm.A[ig,iv] += Float64(-1.00)
        iv,ix_,_ = s2mpj_ii("X11,8",ix_)
        arrset(pb.xnames,iv,"X11,8")
        ig = ig_["N11"]
        pbm.A[ig,iv] += Float64(-1.00)
        ig = ig_["N8"]
        pbm.A[ig,iv] += Float64(1.00)
        iv,ix_,_ = s2mpj_ii("X4,1",ix_)
        arrset(pb.xnames,iv,"X4,1")
        ig = ig_["CNST1"]
        pbm.A[ig,iv] += Float64(20.0)
        ig = ig_["CNST2"]
        pbm.A[ig,iv] += Float64(80.00)
        iv,ix_,_ = s2mpj_ii("X1,4",ix_)
        arrset(pb.xnames,iv,"X1,4")
        ig = ig_["CNST1"]
        pbm.A[ig,iv] += Float64(80.00)
        ig = ig_["CNST2"]
        pbm.A[ig,iv] += Float64(20.00)
        iv,ix_,_ = s2mpj_ii("X4,2",ix_)
        arrset(pb.xnames,iv,"X4,2")
        ig = ig_["CNST3"]
        pbm.A[ig,iv] += Float64(20.0)
        ig = ig_["CNST4"]
        pbm.A[ig,iv] += Float64(80.00)
        iv,ix_,_ = s2mpj_ii("X2,4",ix_)
        arrset(pb.xnames,iv,"X2,4")
        ig = ig_["CNST3"]
        pbm.A[ig,iv] += Float64(80.00)
        ig = ig_["CNST4"]
        pbm.A[ig,iv] += Float64(20.0)
        iv,ix_,_ = s2mpj_ii("X4,3",ix_)
        arrset(pb.xnames,iv,"X4,3")
        ig = ig_["CNST5"]
        pbm.A[ig,iv] += Float64(20.0)
        ig = ig_["CNST6"]
        pbm.A[ig,iv] += Float64(80.00)
        iv,ix_,_ = s2mpj_ii("X3,4",ix_)
        arrset(pb.xnames,iv,"X3,4")
        ig = ig_["CNST5"]
        pbm.A[ig,iv] += Float64(80.00)
        ig = ig_["CNST6"]
        pbm.A[ig,iv] += Float64(20.00)
        iv,ix_,_ = s2mpj_ii("X5,6",ix_)
        arrset(pb.xnames,iv,"X5,6")
        ig = ig_["CNST7"]
        pbm.A[ig,iv] += Float64(20.0)
        ig = ig_["CNST8"]
        pbm.A[ig,iv] += Float64(80.00)
        iv,ix_,_ = s2mpj_ii("X6,5",ix_)
        arrset(pb.xnames,iv,"X6,5")
        ig = ig_["CNST7"]
        pbm.A[ig,iv] += Float64(80.00)
        ig = ig_["CNST8"]
        pbm.A[ig,iv] += Float64(20.0)
        iv,ix_,_ = s2mpj_ii("X5,7",ix_)
        arrset(pb.xnames,iv,"X5,7")
        ig = ig_["CNST9"]
        pbm.A[ig,iv] += Float64(20.0)
        ig = ig_["CNST10"]
        pbm.A[ig,iv] += Float64(80.00)
        iv,ix_,_ = s2mpj_ii("X7,5",ix_)
        arrset(pb.xnames,iv,"X7,5")
        ig = ig_["CNST9"]
        pbm.A[ig,iv] += Float64(80.00)
        ig = ig_["CNST10"]
        pbm.A[ig,iv] += Float64(20.0)
        iv,ix_,_ = s2mpj_ii("X8,9",ix_)
        arrset(pb.xnames,iv,"X8,9")
        ig = ig_["CNST11"]
        pbm.A[ig,iv] += Float64(20.0)
        ig = ig_["CNST12"]
        pbm.A[ig,iv] += Float64(80.00)
        iv,ix_,_ = s2mpj_ii("X9,8",ix_)
        arrset(pb.xnames,iv,"X9,8")
        ig = ig_["CNST11"]
        pbm.A[ig,iv] += Float64(80.00)
        ig = ig_["CNST12"]
        pbm.A[ig,iv] += Float64(20.0)
        iv,ix_,_ = s2mpj_ii("X8,10",ix_)
        arrset(pb.xnames,iv,"X8,10")
        ig = ig_["CNST13"]
        pbm.A[ig,iv] += Float64(20.0)
        ig = ig_["CNST14"]
        pbm.A[ig,iv] += Float64(80.00)
        iv,ix_,_ = s2mpj_ii("X10,8",ix_)
        arrset(pb.xnames,iv,"X10,8")
        ig = ig_["CNST13"]
        pbm.A[ig,iv] += Float64(80.00)
        ig = ig_["CNST14"]
        pbm.A[ig,iv] += Float64(20.0)
        iv,ix_,_ = s2mpj_ii("X8,11",ix_)
        arrset(pb.xnames,iv,"X8,11")
        ig = ig_["CNST15"]
        pbm.A[ig,iv] += Float64(20.0)
        ig = ig_["CNST16"]
        pbm.A[ig,iv] += Float64(80.00)
        iv,ix_,_ = s2mpj_ii("X11,8",ix_)
        arrset(pb.xnames,iv,"X11,8")
        ig = ig_["CNST15"]
        pbm.A[ig,iv] += Float64(80.0)
        ig = ig_["CNST16"]
        pbm.A[ig,iv] += Float64(20.00)
        iv,ix_,_ = s2mpj_ii("X4,5",ix_)
        arrset(pb.xnames,iv,"X4,5")
        ig = ig_["CNST17"]
        pbm.A[ig,iv] += Float64(20.00)
        ig = ig_["CNST18"]
        pbm.A[ig,iv] += Float64(80.0)
        iv,ix_,_ = s2mpj_ii("X5,4",ix_)
        arrset(pb.xnames,iv,"X5,4")
        ig = ig_["CNST17"]
        pbm.A[ig,iv] += Float64(80.0)
        ig = ig_["CNST18"]
        pbm.A[ig,iv] += Float64(20.00)
        iv,ix_,_ = s2mpj_ii("X5,8",ix_)
        arrset(pb.xnames,iv,"X5,8")
        ig = ig_["CNST19"]
        pbm.A[ig,iv] += Float64(20.00)
        ig = ig_["CNST20"]
        pbm.A[ig,iv] += Float64(80.0)
        iv,ix_,_ = s2mpj_ii("X8,5",ix_)
        arrset(pb.xnames,iv,"X8,5")
        ig = ig_["CNST19"]
        pbm.A[ig,iv] += Float64(80.0)
        ig = ig_["CNST20"]
        pbm.A[ig,iv] += Float64(20.00)
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("B"*string(I),ix_)
            arrset(pb.xnames,iv,"B"*string(I))
            ig = ig_["N"*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
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
        pbm.gconst[ig_["N1"]] = Float64(-95.0)
        pbm.gconst[ig_["N2"]] = Float64(-95.0)
        pbm.gconst[ig_["N3"]] = Float64(-19.0)
        pbm.gconst[ig_["N4"]] = Float64(-70.0)
        pbm.gconst[ig_["N5"]] = Float64(-70.0)
        pbm.gconst[ig_["N6"]] = Float64(-19.0)
        pbm.gconst[ig_["N7"]] = Float64(-19.0)
        pbm.gconst[ig_["N8"]] = Float64(-70.0)
        pbm.gconst[ig_["N9"]] = Float64(-19.0)
        pbm.gconst[ig_["N10"]] = Float64(-19.0)
        pbm.gconst[ig_["N11"]] = Float64(-19.0)
        v_["CIJE"] = 999.99
        for I = Int64(v_["1"]):Int64(v_["NLINK-4"])
            pbm.gconst[ig_["CNST"*string(I)]] = Float64(v_["CIJE"])
        end
        v_["CIJE"] = 9999.99
        for I = Int64(v_["NLINK-3"]):Int64(v_["NLINK"])
            pbm.gconst[ig_["CNST"*string(I)]] = Float64(v_["CIJE"])
        end
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xupper[ix_["B1"]] = 99.99
        pb.xupper[ix_["B2"]] = 99.99
        pb.xupper[ix_["B4"]] = 99.99
        pb.xupper[ix_["B5"]] = 99.99
        pb.xupper[ix_["B8"]] = 99.99
        pb.xupper[ix_["B3"]] = 19.99
        pb.xupper[ix_["B6"]] = 19.99
        pb.xupper[ix_["B7"]] = 19.99
        pb.xupper[ix_["B9"]] = 19.99
        pb.xupper[ix_["B10"]] = 19.99
        pb.xupper[ix_["B11"]] = 19.99
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        if haskey(ix_,"X1,4")
            pb.x0[ix_["X1,4"]] = Float64(00.0)
        else
            pb.y0[findfirst(x->x==ig_["X1,4"],pbm.congrps)] = Float64(00.0)
        end
        if haskey(ix_,"X4,1")
            pb.x0[ix_["X4,1"]] = Float64(0.0)
        else
            pb.y0[findfirst(x->x==ig_["X4,1"],pbm.congrps)] = Float64(0.0)
        end
        if haskey(ix_,"B1")
            pb.x0[ix_["B1"]] = Float64(95.0)
        else
            pb.y0[findfirst(x->x==ig_["B1"],pbm.congrps)] = Float64(95.0)
        end
        if haskey(ix_,"X2,4")
            pb.x0[ix_["X2,4"]] = Float64(0.00)
        else
            pb.y0[findfirst(x->x==ig_["X2,4"],pbm.congrps)] = Float64(0.00)
        end
        if haskey(ix_,"X4,2")
            pb.x0[ix_["X4,2"]] = Float64(0.00)
        else
            pb.y0[findfirst(x->x==ig_["X4,2"],pbm.congrps)] = Float64(0.00)
        end
        if haskey(ix_,"B2")
            pb.x0[ix_["B2"]] = Float64(95.0)
        else
            pb.y0[findfirst(x->x==ig_["B2"],pbm.congrps)] = Float64(95.0)
        end
        if haskey(ix_,"X3,4")
            pb.x0[ix_["X3,4"]] = Float64(0.00)
        else
            pb.y0[findfirst(x->x==ig_["X3,4"],pbm.congrps)] = Float64(0.00)
        end
        if haskey(ix_,"X4,3")
            pb.x0[ix_["X4,3"]] = Float64(0.0)
        else
            pb.y0[findfirst(x->x==ig_["X4,3"],pbm.congrps)] = Float64(0.0)
        end
        if haskey(ix_,"B3")
            pb.x0[ix_["B3"]] = Float64(19.0)
        else
            pb.y0[findfirst(x->x==ig_["B3"],pbm.congrps)] = Float64(19.0)
        end
        if haskey(ix_,"B4")
            pb.x0[ix_["B4"]] = Float64(70.0)
        else
            pb.y0[findfirst(x->x==ig_["B4"],pbm.congrps)] = Float64(70.0)
        end
        if haskey(ix_,"X5,4")
            pb.x0[ix_["X5,4"]] = Float64(0.0)
        else
            pb.y0[findfirst(x->x==ig_["X5,4"],pbm.congrps)] = Float64(0.0)
        end
        if haskey(ix_,"X4,5")
            pb.x0[ix_["X4,5"]] = Float64(00.0)
        else
            pb.y0[findfirst(x->x==ig_["X4,5"],pbm.congrps)] = Float64(00.0)
        end
        if haskey(ix_,"X6,5")
            pb.x0[ix_["X6,5"]] = Float64(0.0)
        else
            pb.y0[findfirst(x->x==ig_["X6,5"],pbm.congrps)] = Float64(0.0)
        end
        if haskey(ix_,"X5,6")
            pb.x0[ix_["X5,6"]] = Float64(0.0)
        else
            pb.y0[findfirst(x->x==ig_["X5,6"],pbm.congrps)] = Float64(0.0)
        end
        if haskey(ix_,"X7,5")
            pb.x0[ix_["X7,5"]] = Float64(0.0)
        else
            pb.y0[findfirst(x->x==ig_["X7,5"],pbm.congrps)] = Float64(0.0)
        end
        if haskey(ix_,"X5,7")
            pb.x0[ix_["X5,7"]] = Float64(0.0)
        else
            pb.y0[findfirst(x->x==ig_["X5,7"],pbm.congrps)] = Float64(0.0)
        end
        if haskey(ix_,"B5")
            pb.x0[ix_["B5"]] = Float64(70.0)
        else
            pb.y0[findfirst(x->x==ig_["B5"],pbm.congrps)] = Float64(70.0)
        end
        if haskey(ix_,"B6")
            pb.x0[ix_["B6"]] = Float64(19.0)
        else
            pb.y0[findfirst(x->x==ig_["B6"],pbm.congrps)] = Float64(19.0)
        end
        if haskey(ix_,"B7")
            pb.x0[ix_["B7"]] = Float64(19.0)
        else
            pb.y0[findfirst(x->x==ig_["B7"],pbm.congrps)] = Float64(19.0)
        end
        if haskey(ix_,"X8,5")
            pb.x0[ix_["X8,5"]] = Float64(0.0)
        else
            pb.y0[findfirst(x->x==ig_["X8,5"],pbm.congrps)] = Float64(0.0)
        end
        if haskey(ix_,"X5,8")
            pb.x0[ix_["X5,8"]] = Float64(0.0)
        else
            pb.y0[findfirst(x->x==ig_["X5,8"],pbm.congrps)] = Float64(0.0)
        end
        if haskey(ix_,"X9,8")
            pb.x0[ix_["X9,8"]] = Float64(0.0)
        else
            pb.y0[findfirst(x->x==ig_["X9,8"],pbm.congrps)] = Float64(0.0)
        end
        if haskey(ix_,"X8,9")
            pb.x0[ix_["X8,9"]] = Float64(0.0)
        else
            pb.y0[findfirst(x->x==ig_["X8,9"],pbm.congrps)] = Float64(0.0)
        end
        if haskey(ix_,"X10,8")
            pb.x0[ix_["X10,8"]] = Float64(0.0)
        else
            pb.y0[findfirst(x->x==ig_["X10,8"],pbm.congrps)] = Float64(0.0)
        end
        if haskey(ix_,"X8,10")
            pb.x0[ix_["X8,10"]] = Float64(0.0)
        else
            pb.y0[findfirst(x->x==ig_["X8,10"],pbm.congrps)] = Float64(0.0)
        end
        if haskey(ix_,"X11,8")
            pb.x0[ix_["X11,8"]] = Float64(0.0)
        else
            pb.y0[findfirst(x->x==ig_["X11,8"],pbm.congrps)] = Float64(0.0)
        end
        if haskey(ix_,"X8,11")
            pb.x0[ix_["X8,11"]] = Float64(0.0)
        else
            pb.y0[findfirst(x->x==ig_["X8,11"],pbm.congrps)] = Float64(0.0)
        end
        if haskey(ix_,"B8")
            pb.x0[ix_["B8"]] = Float64(70.0)
        else
            pb.y0[findfirst(x->x==ig_["B8"],pbm.congrps)] = Float64(70.0)
        end
        if haskey(ix_,"B9")
            pb.x0[ix_["B9"]] = Float64(19.0)
        else
            pb.y0[findfirst(x->x==ig_["B9"],pbm.congrps)] = Float64(19.0)
        end
        if haskey(ix_,"B10")
            pb.x0[ix_["B10"]] = Float64(19.0)
        else
            pb.y0[findfirst(x->x==ig_["B10"],pbm.congrps)] = Float64(19.0)
        end
        if haskey(ix_,"B11")
            pb.x0[ix_["B11"]] = Float64(19.0)
        else
            pb.y0[findfirst(x->x==ig_["B11"],pbm.congrps)] = Float64(19.0)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eBETA1", iet_)
        loaset(elftv,it,1,"V")
        it,iet_,_ = s2mpj_ii( "eBETA2", iet_)
        loaset(elftv,it,1,"V")
        it,iet_,_ = s2mpj_ii( "eCOMA1", iet_)
        loaset(elftv,it,1,"V")
        loaset(elftv,it,2,"W")
        it,iet_,_ = s2mpj_ii( "eCOMA2", iet_)
        loaset(elftv,it,1,"V")
        loaset(elftv,it,2,"W")
        it,iet_,_ = s2mpj_ii( "eCOMB1", iet_)
        loaset(elftv,it,1,"V")
        loaset(elftv,it,2,"W")
        it,iet_,_ = s2mpj_ii( "eCOMB2", iet_)
        loaset(elftv,it,1,"V")
        loaset(elftv,it,2,"W")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "EB1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eBETA1")
        arrset(ielftype,ie,iet_["eBETA1"])
        vname = "B1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EB2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eBETA1")
        arrset(ielftype,ie,iet_["eBETA1"])
        vname = "B2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EB3"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eBETA2")
        arrset(ielftype,ie,iet_["eBETA2"])
        vname = "B3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EB4"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eBETA1")
        arrset(ielftype,ie,iet_["eBETA1"])
        vname = "B4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EB5"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eBETA1")
        arrset(ielftype,ie,iet_["eBETA1"])
        vname = "B5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EB6"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eBETA2")
        arrset(ielftype,ie,iet_["eBETA2"])
        vname = "B6"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EB7"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eBETA2")
        arrset(ielftype,ie,iet_["eBETA2"])
        vname = "B7"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EB8"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eBETA1")
        arrset(ielftype,ie,iet_["eBETA1"])
        vname = "B8"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EB9"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eBETA2")
        arrset(ielftype,ie,iet_["eBETA2"])
        vname = "B9"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EB10"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eBETA2")
        arrset(ielftype,ie,iet_["eBETA2"])
        vname = "B10"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EB11"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eBETA2")
        arrset(ielftype,ie,iet_["eBETA2"])
        vname = "B11"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        for I = Int64(v_["1"]):Int64(v_["P1"])
            ename = "EGA"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eCOMA1")
            arrset(ielftype,ie,iet_["eCOMA1"])
            vname = "X"*string(Int64(v_["4C"]))*","*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(I)*","*string(Int64(v_["4C"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="W",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "EGB"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eCOMB1")
            arrset(ielftype,ie,iet_["eCOMB1"])
            vname = "X"*string(Int64(v_["4C"]))*","*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(I)*","*string(Int64(v_["4C"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="W",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            v_["I+3"] = 3+I
            ename = "EGA"*string(Int64(v_["I+3"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eCOMA1")
            arrset(ielftype,ie,iet_["eCOMA1"])
            ename = "EGA"*string(Int64(v_["I+3"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X"*string(I)*","*string(Int64(v_["4C"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "EGA"*string(Int64(v_["I+3"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X"*string(Int64(v_["4C"]))*","*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="W",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "EGB"*string(Int64(v_["I+3"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eCOMB1")
            arrset(ielftype,ie,iet_["eCOMB1"])
            ename = "EGB"*string(Int64(v_["I+3"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X"*string(I)*","*string(Int64(v_["4C"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "EGB"*string(Int64(v_["I+3"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X"*string(Int64(v_["4C"]))*","*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="W",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            v_["I+6"] = 6+I
            v_["I+8"] = 8+I
            ename = "EGA"*string(Int64(v_["I+6"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eCOMA1")
            arrset(ielftype,ie,iet_["eCOMA1"])
            ename = "EGA"*string(Int64(v_["I+6"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X"*string(Int64(v_["8C"]))*","*string(Int64(v_["I+8"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "EGA"*string(Int64(v_["I+6"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X"*string(Int64(v_["I+8"]))*","*string(Int64(v_["8C"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="W",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "EGB"*string(Int64(v_["I+6"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eCOMB1")
            arrset(ielftype,ie,iet_["eCOMB1"])
            ename = "EGB"*string(Int64(v_["I+6"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X"*string(Int64(v_["8C"]))*","*string(Int64(v_["I+8"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "EGB"*string(Int64(v_["I+6"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X"*string(Int64(v_["I+8"]))*","*string(Int64(v_["8C"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="W",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            v_["I+9"] = 9+I
            ename = "EGA"*string(Int64(v_["I+9"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eCOMA1")
            arrset(ielftype,ie,iet_["eCOMA1"])
            ename = "EGA"*string(Int64(v_["I+9"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X"*string(Int64(v_["I+8"]))*","*string(Int64(v_["8C"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "EGA"*string(Int64(v_["I+9"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X"*string(Int64(v_["8C"]))*","*string(Int64(v_["I+8"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="W",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "EGB"*string(Int64(v_["I+9"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eCOMB1")
            arrset(ielftype,ie,iet_["eCOMB1"])
            ename = "EGB"*string(Int64(v_["I+9"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X"*string(Int64(v_["I+8"]))*","*string(Int64(v_["8C"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "EGB"*string(Int64(v_["I+9"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X"*string(Int64(v_["8C"]))*","*string(Int64(v_["I+8"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="W",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        for I = Int64(v_["6C"]):Int64(v_["7C"])
            v_["I2"] = 2*I
            v_["I2+1"] = 1+v_["I2"]
            ename = "EGA"*string(Int64(v_["I2+1"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eCOMA1")
            arrset(ielftype,ie,iet_["eCOMA1"])
            ename = "EGA"*string(Int64(v_["I2+1"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X"*string(Int64(v_["5C"]))*","*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "EGA"*string(Int64(v_["I2+1"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X"*string(I)*","*string(Int64(v_["5C"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="W",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "EGB"*string(Int64(v_["I2+1"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eCOMB1")
            arrset(ielftype,ie,iet_["eCOMB1"])
            ename = "EGB"*string(Int64(v_["I2+1"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X"*string(Int64(v_["5C"]))*","*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "EGB"*string(Int64(v_["I2+1"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X"*string(I)*","*string(Int64(v_["5C"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="W",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            v_["I2+2"] = 2+v_["I2"]
            ename = "EGA"*string(Int64(v_["I2+2"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eCOMA1")
            arrset(ielftype,ie,iet_["eCOMA1"])
            ename = "EGA"*string(Int64(v_["I2+2"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X"*string(I)*","*string(Int64(v_["5C"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "EGA"*string(Int64(v_["I2+2"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X"*string(Int64(v_["5C"]))*","*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="W",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "EGB"*string(Int64(v_["I2+2"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eCOMB1")
            arrset(ielftype,ie,iet_["eCOMB1"])
            ename = "EGB"*string(Int64(v_["I2+2"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X"*string(I)*","*string(Int64(v_["5C"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "EGB"*string(Int64(v_["I2+2"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X"*string(Int64(v_["5C"]))*","*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="W",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        ename = "EGA17"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eCOMA2")
        arrset(ielftype,ie,iet_["eCOMA2"])
        vname = "X5,4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X4,5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="W",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EGB17"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eCOMB2")
        arrset(ielftype,ie,iet_["eCOMB2"])
        vname = "X5,4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X4,5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="W",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EGA18"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eCOMA2")
        arrset(ielftype,ie,iet_["eCOMA2"])
        vname = "X4,5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X5,4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="W",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EGB18"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eCOMB2")
        arrset(ielftype,ie,iet_["eCOMB2"])
        vname = "X4,5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X5,4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="W",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EGA19"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eCOMA2")
        arrset(ielftype,ie,iet_["eCOMA2"])
        vname = "X5,8"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X8,5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="W",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EGB19"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eCOMB2")
        arrset(ielftype,ie,iet_["eCOMB2"])
        vname = "X5,8"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X8,5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="W",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EGA20"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eCOMA2")
        arrset(ielftype,ie,iet_["eCOMA2"])
        vname = "X8,5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X5,8"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="W",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EGB20"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eCOMB2")
        arrset(ielftype,ie,iet_["eCOMB2"])
        vname = "X8,5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X5,8"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="W",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            ig = ig_["F"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["EB"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        for I = Int64(v_["1"]):Int64(v_["NLINK"])
            ig = ig_["GB"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["EGB"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            ig = ig_["GA"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["EGA"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-OLR2-MN-31-31"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm

# **********************
#  SET UP THE ELEMENTS *
#  ROUTINE             *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "e_globs"

        pbm = args[1]
        arrset(pbm.efpar,1,80.0)
        arrset(pbm.efpar,2,20.0)
        return pbm

    elseif action == "eBETA2"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        CB = 20.0
        f_   = EV_[1]/(CB-EV_[1])
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = CB/((CB-EV_[1])^2)
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 2*CB/((CB-EV_[1])^3)
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eBETA1"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        CB = 100.0
        f_   = EV_[1]/(CB-EV_[1])
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = CB/((CB-EV_[1])^2)
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 2*CB/((CB-EV_[1])^3)
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eCOMA1"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        CIJ = 1000.0
        f_   = EV_[1]/(CIJ-(pbm.efpar[1]*EV_[1]+pbm.efpar[2]*EV_[2]))
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1]  = (
                  (CIJ-pbm.efpar[2]*EV_[2])/(CIJ-(pbm.efpar[1]*EV_[1]+pbm.efpar[2]*EV_[2]))^2)
            g_[2]  = (
                  EV_[1]*pbm.efpar[2]/(CIJ-(pbm.efpar[1]*EV_[1]+pbm.efpar[2]*EV_[2]))^2)
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1]  = (
                      2*pbm.efpar[1]*(CIJ-pbm.efpar[2]*EV_[2])/(CIJ-(pbm.efpar[1]*EV_[1]+pbm.efpar[2]*EV_[2]))^3)
                H_[1,2] = (pbm.efpar[2]*(CIJ+pbm.efpar[1]*EV_[1]-pbm.efpar[2]*EV_[2])/
                     (CIJ-(pbm.efpar[1]*EV_[1]+pbm.efpar[2]*EV_[2]))^3)
                H_[2,1] = H_[1,2]
                H_[2,2]  = (
                      2*pbm.efpar[2]*pbm.efpar[2]*EV_[1]/(CIJ-(pbm.efpar[1]*EV_[1]+pbm.efpar[2]*EV_[2]))^3)
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eCOMB1"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        CIJ = 1000.0
        f_   = EV_[1]/(CIJ-(pbm.efpar[2]*EV_[1]+pbm.efpar[1]*EV_[2]))
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1]  = (
                  (CIJ-pbm.efpar[1]*EV_[2])/(CIJ-(pbm.efpar[2]*EV_[1]+pbm.efpar[1]*EV_[2]))^2)
            g_[2]  = (
                  EV_[1]*pbm.efpar[1]/(CIJ-(pbm.efpar[2]*EV_[1]+pbm.efpar[1]*EV_[2]))^2)
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1]  = (
                      2*pbm.efpar[2]*(CIJ-pbm.efpar[1]*EV_[2])/(CIJ-(pbm.efpar[2]*EV_[1]+pbm.efpar[1]*EV_[2]))^3)
                H_[1,2] = (pbm.efpar[1]*(CIJ+pbm.efpar[2]*EV_[1]-pbm.efpar[1]*EV_[2])/
                     (CIJ-(pbm.efpar[2]*EV_[1]+pbm.efpar[1]*EV_[2]))^3)
                H_[2,1] = H_[1,2]
                H_[2,2]  = (
                      2*pbm.efpar[1]*pbm.efpar[1]*EV_[1]/(CIJ-(pbm.efpar[2]*EV_[1]+pbm.efpar[1]*EV_[2]))^3)
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eCOMA2"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        CIJ = 10000.0
        f_   = EV_[1]/(CIJ-(pbm.efpar[1]*EV_[1]+pbm.efpar[2]*EV_[2]))
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1]  = (
                  (CIJ-pbm.efpar[2]*EV_[2])/(CIJ-(pbm.efpar[1]*EV_[1]+pbm.efpar[2]*EV_[2]))^2)
            g_[2]  = (
                  EV_[1]*pbm.efpar[2]/(CIJ-(pbm.efpar[1]*EV_[1]+pbm.efpar[2]*EV_[2]))^2)
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1]  = (
                      2*pbm.efpar[1]*(CIJ-pbm.efpar[2]*EV_[2])/(CIJ-(pbm.efpar[1]*EV_[1]+pbm.efpar[2]*EV_[2]))^3)
                H_[1,2] = (pbm.efpar[2]*(CIJ+pbm.efpar[1]*EV_[1]-pbm.efpar[2]*EV_[2])/
                     (CIJ-(pbm.efpar[1]*EV_[1]+pbm.efpar[2]*EV_[2]))^3)
                H_[2,1] = H_[1,2]
                H_[2,2]  = (
                      2*pbm.efpar[2]*pbm.efpar[2]*EV_[1]/(CIJ-(pbm.efpar[1]*EV_[1]+pbm.efpar[2]*EV_[2]))^3)
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eCOMB2"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        CIJ = 10000.0
        f_   = EV_[1]/(CIJ-(pbm.efpar[2]*EV_[1]+pbm.efpar[1]*EV_[2]))
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1]  = (
                  (CIJ-pbm.efpar[1]*EV_[2])/(CIJ-(pbm.efpar[2]*EV_[1]+pbm.efpar[1]*EV_[2]))^2)
            g_[2]  = (
                  EV_[1]*pbm.efpar[1]/(CIJ-(pbm.efpar[2]*EV_[1]+pbm.efpar[1]*EV_[2]))^2)
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1]  = (
                      2*pbm.efpar[2]*(CIJ-pbm.efpar[1]*EV_[2])/(CIJ-(pbm.efpar[2]*EV_[1]+pbm.efpar[1]*EV_[2]))^3)
                H_[1,2] = (pbm.efpar[1]*(CIJ+pbm.efpar[2]*EV_[1]-pbm.efpar[1]*EV_[2])/
                     (CIJ-(pbm.efpar[2]*EV_[1]+pbm.efpar[1]*EV_[2]))^3)
                H_[2,1] = H_[1,2]
                H_[2,2]  = (
                      2*pbm.efpar[1]*pbm.efpar[1]*EV_[1]/(CIJ-(pbm.efpar[2]*EV_[1]+pbm.efpar[1]*EV_[2]))^3)
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
            pbm.has_globs = [2,0]
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

