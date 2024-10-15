function LAUNCH(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    The objective function to be minimized represents the total cost of
#    the development and launching of a 3 stages space launching vehicle.
#    Constraints are imposed on physical interrelations between the variables
#    and performance.
# 
#    The problem is highly non-convex. 
# 
#    Source:
#    B. Rush, J. Bracken and G. McCormick,
#    "A nonliner programming model for launch vehicle design and costing",
#    Operations Research, pp. 185-210, 1967.
# 
#    SIF input: P. Driscoll, Virginia Tech., April 1993,
#               corrected and simplified by Ph. L. Toint, May 1993.
# 
#    classification = "C-OOR2-MY-25-28"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "LAUNCH"

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
        iv,ix_,_ = s2mpj_ii("AW1",ix_)
        arrset(pb.xnames,iv,"AW1")
        iv,ix_,_ = s2mpj_ii("IW1",ix_)
        arrset(pb.xnames,iv,"IW1")
        iv,ix_,_ = s2mpj_ii("MF1",ix_)
        arrset(pb.xnames,iv,"MF1")
        iv,ix_,_ = s2mpj_ii("TT1",ix_)
        arrset(pb.xnames,iv,"TT1")
        iv,ix_,_ = s2mpj_ii("PW1",ix_)
        arrset(pb.xnames,iv,"PW1")
        iv,ix_,_ = s2mpj_ii("ET1",ix_)
        arrset(pb.xnames,iv,"ET1")
        iv,ix_,_ = s2mpj_ii("S1L",ix_)
        arrset(pb.xnames,iv,"S1L")
        iv,ix_,_ = s2mpj_ii("AW2",ix_)
        arrset(pb.xnames,iv,"AW2")
        iv,ix_,_ = s2mpj_ii("IW2",ix_)
        arrset(pb.xnames,iv,"IW2")
        iv,ix_,_ = s2mpj_ii("MF2",ix_)
        arrset(pb.xnames,iv,"MF2")
        iv,ix_,_ = s2mpj_ii("TT2",ix_)
        arrset(pb.xnames,iv,"TT2")
        iv,ix_,_ = s2mpj_ii("PW2",ix_)
        arrset(pb.xnames,iv,"PW2")
        iv,ix_,_ = s2mpj_ii("ET2",ix_)
        arrset(pb.xnames,iv,"ET2")
        iv,ix_,_ = s2mpj_ii("S2L",ix_)
        arrset(pb.xnames,iv,"S2L")
        iv,ix_,_ = s2mpj_ii("AW3",ix_)
        arrset(pb.xnames,iv,"AW3")
        iv,ix_,_ = s2mpj_ii("IW3",ix_)
        arrset(pb.xnames,iv,"IW3")
        iv,ix_,_ = s2mpj_ii("MF3",ix_)
        arrset(pb.xnames,iv,"MF3")
        iv,ix_,_ = s2mpj_ii("TT3",ix_)
        arrset(pb.xnames,iv,"TT3")
        iv,ix_,_ = s2mpj_ii("PW3",ix_)
        arrset(pb.xnames,iv,"PW3")
        iv,ix_,_ = s2mpj_ii("ET3",ix_)
        arrset(pb.xnames,iv,"ET3")
        iv,ix_,_ = s2mpj_ii("S3L",ix_)
        arrset(pb.xnames,iv,"S3L")
        iv,ix_,_ = s2mpj_ii("INW",ix_)
        arrset(pb.xnames,iv,"INW")
        iv,ix_,_ = s2mpj_ii("BT1",ix_)
        arrset(pb.xnames,iv,"BT1")
        iv,ix_,_ = s2mpj_ii("BT2",ix_)
        arrset(pb.xnames,iv,"BT2")
        iv,ix_,_ = s2mpj_ii("BT3",ix_)
        arrset(pb.xnames,iv,"BT3")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("STA1",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["ET1"]
        pbm.A[ig,iv] += Float64(0.0002587)
        arrset(pbm.gscale,ig,Float64(1.0e+8))
        ig,ig_,_ = s2mpj_ii("STA2",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["ET2"]
        pbm.A[ig,iv] += Float64(0.0002587)
        arrset(pbm.gscale,ig,Float64(1.0e+8))
        ig,ig_,_ = s2mpj_ii("STA3",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["ET3"]
        pbm.A[ig,iv] += Float64(0.001958)
        arrset(pbm.gscale,ig,Float64(1.0e+8))
        ig,ig_,_ = s2mpj_ii("INST",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["INW"]
        pbm.A[ig,iv] += Float64(47.040096)
        arrset(pbm.gscale,ig,Float64(1.0e+8))
        ig,ig_,_ = s2mpj_ii("LAUN",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["PW1"]
        pbm.A[ig,iv] += Float64(0.003)
        iv = ix_["PW2"]
        pbm.A[ig,iv] += Float64(0.003)
        iv = ix_["PW3"]
        pbm.A[ig,iv] += Float64(0.003)
        arrset(pbm.gscale,ig,Float64(39215686.0))
        ig,ig_,_ = s2mpj_ii("SGTH1",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"SGTH1")
        iv = ix_["AW1"]
        pbm.A[ig,iv] += Float64(2.0)
        iv = ix_["IW1"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("SGTH3",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"SGTH3")
        iv = ix_["IW2"]
        pbm.A[ig,iv] += Float64(0.6)
        iv = ix_["AW2"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("SGTH5",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"SGTH5")
        iv = ix_["IW3"]
        pbm.A[ig,iv] += Float64(0.7)
        iv = ix_["AW3"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("SGTH2",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"SGTH2")
        iv = ix_["ET1"]
        pbm.A[ig,iv] += Float64(5.0)
        iv = ix_["TT1"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("SGTH4",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"SGTH4")
        iv = ix_["ET2"]
        pbm.A[ig,iv] += Float64(5.0)
        iv = ix_["TT2"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("SGTH6",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"SGTH6")
        iv = ix_["TT3"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["ET3"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("SGSI1A",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"SGSI1A")
        iv = ix_["PW1"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["IW1"]
        pbm.A[ig,iv] += Float64(-12.0)
        ig,ig_,_ = s2mpj_ii("SGSI1B",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"SGSI1B")
        iv = ix_["PW1"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["IW1"]
        pbm.A[ig,iv] += Float64(-17.0)
        ig,ig_,_ = s2mpj_ii("SGSI2A",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"SGSI2A")
        iv = ix_["PW2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["IW2"]
        pbm.A[ig,iv] += Float64(-10.0)
        ig,ig_,_ = s2mpj_ii("SGSI2B",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"SGSI2B")
        iv = ix_["PW2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["IW2"]
        pbm.A[ig,iv] += Float64(-13.0)
        ig,ig_,_ = s2mpj_ii("SGSI3A",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"SGSI3A")
        iv = ix_["PW3"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["IW3"]
        pbm.A[ig,iv] += Float64(-7.0)
        ig,ig_,_ = s2mpj_ii("SGSI3B",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"SGSI3B")
        iv = ix_["PW3"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["IW3"]
        pbm.A[ig,iv] += Float64(-10.0)
        ig,ig_,_ = s2mpj_ii("TTIW1A",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"TTIW1A")
        iv = ix_["TT1"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["IW1"]
        pbm.A[ig,iv] += Float64(-1.2)
        iv = ix_["PW1"]
        pbm.A[ig,iv] += Float64(-1.2)
        iv = ix_["IW2"]
        pbm.A[ig,iv] += Float64(-1.2)
        iv = ix_["PW2"]
        pbm.A[ig,iv] += Float64(-1.2)
        iv = ix_["IW3"]
        pbm.A[ig,iv] += Float64(-1.2)
        iv = ix_["PW3"]
        pbm.A[ig,iv] += Float64(-1.2)
        iv = ix_["INW"]
        pbm.A[ig,iv] += Float64(-1.2)
        ig,ig_,_ = s2mpj_ii("TTIW1B",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"TTIW1B")
        iv = ix_["TT1"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["IW1"]
        pbm.A[ig,iv] += Float64(-1.4)
        iv = ix_["PW1"]
        pbm.A[ig,iv] += Float64(-1.4)
        iv = ix_["IW2"]
        pbm.A[ig,iv] += Float64(-1.4)
        iv = ix_["PW2"]
        pbm.A[ig,iv] += Float64(-1.4)
        iv = ix_["IW3"]
        pbm.A[ig,iv] += Float64(-1.4)
        iv = ix_["PW3"]
        pbm.A[ig,iv] += Float64(-1.4)
        iv = ix_["INW"]
        pbm.A[ig,iv] += Float64(-1.4)
        ig,ig_,_ = s2mpj_ii("TTIW2A",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"TTIW2A")
        iv = ix_["TT2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["IW2"]
        pbm.A[ig,iv] += Float64(-0.6)
        iv = ix_["PW2"]
        pbm.A[ig,iv] += Float64(-0.6)
        iv = ix_["IW3"]
        pbm.A[ig,iv] += Float64(-0.6)
        iv = ix_["PW3"]
        pbm.A[ig,iv] += Float64(-0.6)
        iv = ix_["INW"]
        pbm.A[ig,iv] += Float64(-0.6)
        ig,ig_,_ = s2mpj_ii("TTIW2B",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"TTIW2B")
        iv = ix_["TT2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["IW2"]
        pbm.A[ig,iv] += Float64(-0.75)
        iv = ix_["PW2"]
        pbm.A[ig,iv] += Float64(-0.75)
        iv = ix_["IW3"]
        pbm.A[ig,iv] += Float64(-0.75)
        iv = ix_["PW3"]
        pbm.A[ig,iv] += Float64(-0.75)
        iv = ix_["INW"]
        pbm.A[ig,iv] += Float64(-0.75)
        ig,ig_,_ = s2mpj_ii("TTIW3A",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"TTIW3A")
        iv = ix_["TT3"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["IW3"]
        pbm.A[ig,iv] += Float64(-0.7)
        iv = ix_["PW3"]
        pbm.A[ig,iv] += Float64(-0.7)
        iv = ix_["INW"]
        pbm.A[ig,iv] += Float64(-0.7)
        ig,ig_,_ = s2mpj_ii("TTIW3B",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"TTIW3B")
        iv = ix_["TT3"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["IW3"]
        pbm.A[ig,iv] += Float64(-0.9)
        iv = ix_["PW3"]
        pbm.A[ig,iv] += Float64(-0.9)
        iv = ix_["INW"]
        pbm.A[ig,iv] += Float64(-0.9)
        ig,ig_,_ = s2mpj_ii("SMF1",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"SMF1")
        iv = ix_["MF1"]
        pbm.A[ig,iv] += Float64(20.0)
        iv = ix_["IW1"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["IW2"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["PW2"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["IW3"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["PW3"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["INW"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("SMF2",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"SMF2")
        iv = ix_["MF2"]
        pbm.A[ig,iv] += Float64(20.0)
        iv = ix_["IW2"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["IW3"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["PW3"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["INW"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("SMF3",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"SMF3")
        iv = ix_["MF3"]
        pbm.A[ig,iv] += Float64(20.0)
        iv = ix_["IW3"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["INW"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("SI1A",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"SI1A")
        iv = ix_["PW1"]
        pbm.A[ig,iv] += Float64(-240.0)
        ig,ig_,_ = s2mpj_ii("SI1B",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"SI1B")
        iv = ix_["PW1"]
        pbm.A[ig,iv] += Float64(-290.0)
        ig,ig_,_ = s2mpj_ii("SI2A",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"SI2A")
        iv = ix_["PW2"]
        pbm.A[ig,iv] += Float64(-240.0)
        ig,ig_,_ = s2mpj_ii("SI2B",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"SI2B")
        iv = ix_["PW2"]
        pbm.A[ig,iv] += Float64(-290.0)
        ig,ig_,_ = s2mpj_ii("SI3A",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"SI3A")
        iv = ix_["PW3"]
        pbm.A[ig,iv] += Float64(-340.0)
        ig,ig_,_ = s2mpj_ii("SI3B",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"SI3B")
        iv = ix_["PW3"]
        pbm.A[ig,iv] += Float64(-375.0)
        ig,ig_,_ = s2mpj_ii("GLGCON",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"GLGCON")
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
        pbm.gconst[ig_["STA1"]] = Float64(247.963)
        pbm.gconst[ig_["STA2"]] = Float64(247.963)
        pbm.gconst[ig_["STA3"]] = Float64(32.591)
        pbm.gconst[ig_["INST"]] = Float64(35.5)
        pbm.gconst[ig_["TTIW1A"]] = Float64(24.0)
        pbm.gconst[ig_["TTIW1B"]] = Float64(28.0)
        pbm.gconst[ig_["TTIW2A"]] = Float64(12.0)
        pbm.gconst[ig_["TTIW2B"]] = Float64(15.0)
        pbm.gconst[ig_["TTIW3A"]] = Float64(14.0)
        pbm.gconst[ig_["TTIW3B"]] = Float64(18.0)
        pbm.gconst[ig_["SMF1"]] = Float64(20.0)
        pbm.gconst[ig_["SMF2"]] = Float64(20.0)
        pbm.gconst[ig_["SMF3"]] = Float64(20.0)
        pbm.gconst[ig_["GLGCON"]] = Float64(35000.0)
        #%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange = Vector{Float64}(undef,ngrp)
        grange[legrps,1] = fill(Inf,pb.nle)
        grange[gegrps,1] = fill(Inf,pb.nge)
        arrset(grange,ig_["GLGCON"],Float64(15000.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = fill(1.0e-8,pb.n)
        pb.xupper = fill(1.0e+4,pb.n)
        pb.xlower[ix_["S1L"]] = 125.0
        pb.xupper[ix_["S1L"]] = 150.0
        pb.xlower[ix_["S2L"]] = 75.0
        pb.xupper[ix_["S2L"]] = 100.0
        pb.xlower[ix_["S3L"]] = 50.0
        pb.xupper[ix_["S3L"]] = 70.0
        pb.xlower[ix_["MF1"]] = 0.25
        pb.xupper[ix_["MF1"]] = 0.30
        pb.xlower[ix_["MF2"]] = 0.24
        pb.xupper[ix_["MF2"]] = 0.29
        pb.xlower[ix_["MF3"]] = 0.16
        pb.xupper[ix_["MF3"]] = 0.21
        pb.xlower[ix_["INW"]] = 2.5
        pb.xupper[ix_["INW"]] = 4.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        pb.x0[ix_["AW1"]] = Float64(68.0)
        pb.x0[ix_["IW1"]] = Float64(136.0)
        pb.x0[ix_["MF1"]] = Float64(0.29988744)
        pb.x0[ix_["TT1"]] = Float64(3733.0)
        pb.x0[ix_["PW1"]] = Float64(2177.0)
        pb.x0[ix_["ET1"]] = Float64(746.6)
        pb.x0[ix_["S1L"]] = Float64(125.0)
        pb.x0[ix_["AW2"]] = Float64(28.2)
        pb.x0[ix_["IW2"]] = Float64(47.0)
        pb.x0[ix_["MF2"]] = Float64(0.28939109)
        pb.x0[ix_["TT2"]] = Float64(480.0)
        pb.x0[ix_["PW2"]] = Float64(566.0)
        pb.x0[ix_["ET2"]] = Float64(96.0)
        pb.x0[ix_["S2L"]] = Float64(75.0)
        pb.x0[ix_["AW3"]] = Float64(11.2)
        pb.x0[ix_["IW3"]] = Float64(16.0)
        pb.x0[ix_["MF3"]] = Float64(0.20980926)
        pb.x0[ix_["TT3"]] = Float64(129.0)
        pb.x0[ix_["PW3"]] = Float64(145.0)
        pb.x0[ix_["ET3"]] = Float64(129.0)
        pb.x0[ix_["S3L"]] = Float64(50.0)
        pb.x0[ix_["INW"]] = Float64(2.5)
        pb.x0[ix_["BT1"]] = Float64(155.0)
        pb.x0[ix_["BT2"]] = Float64(314.0)
        pb.x0[ix_["BT3"]] = Float64(403.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "ePROD1", iet_)
        loaset(elftv,it,1,"VA")
        loaset(elftv,it,2,"VB")
        loaset(elftv,it,3,"VC")
        loaset(elftv,it,4,"VD")
        loaset(elftv,it,5,"VE")
        it,iet_,_ = s2mpj_ii( "ePOWER", iet_)
        loaset(elftv,it,1,"XX")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"PWR")
        loaset(elftp,it,2,"SC")
        it,iet_,_ = s2mpj_ii( "ePROD2", iet_)
        loaset(elftv,it,1,"VA")
        loaset(elftv,it,2,"VB")
        loaset(elftv,it,3,"VC")
        loaset(elftv,it,4,"VD")
        it,iet_,_ = s2mpj_ii( "eX7Y", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y1")
        loaset(elftv,it,3,"Y2")
        loaset(elftv,it,4,"Y3")
        loaset(elftv,it,5,"Y4")
        loaset(elftv,it,6,"Y5")
        loaset(elftv,it,7,"Y6")
        loaset(elftv,it,8,"Y7")
        it,iet_,_ = s2mpj_ii( "eX5Y", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y1")
        loaset(elftv,it,3,"Y2")
        loaset(elftv,it,4,"Y3")
        loaset(elftv,it,5,"Y4")
        loaset(elftv,it,6,"Y5")
        it,iet_,_ = s2mpj_ii( "eX3Y", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y1")
        loaset(elftv,it,3,"Y2")
        loaset(elftv,it,4,"Y3")
        it,iet_,_ = s2mpj_ii( "en2PR", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        it,iet_,_ = s2mpj_ii( "eBIG1", iet_)
        loaset(elftv,it,1,"LH")
        loaset(elftv,it,2,"TH")
        loaset(elftv,it,3,"LL")
        loaset(elftv,it,4,"V1")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "XPROD1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD1")
        arrset(ielftype,ie,iet_["ePROD1"])
        vname = "AW1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="VA",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "IW1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="VB",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "MF1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="VC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TT1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="VD",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PW1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="VE",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "XPF"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePOWER")
        arrset(ielftype,ie,iet_["ePOWER"])
        vname = "ET1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="PWR",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-0.146))
        posep = findfirst(x->x=="SC",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1000.0))
        ename = "XPG"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePOWER")
        arrset(ielftype,ie,iet_["ePOWER"])
        vname = "ET1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="PWR",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.648))
        posep = findfirst(x->x=="SC",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1000.0))
        ename = "XPROD2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD2")
        arrset(ielftype,ie,iet_["ePROD2"])
        vname = "AW1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="VA",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "MF1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="VB",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PW1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="VC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "S1L"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="VD",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "XPL"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePOWER")
        arrset(ielftype,ie,iet_["ePOWER"])
        vname = "ET1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="PWR",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.736))
        posep = findfirst(x->x=="SC",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1000.0))
        ename = "XPM"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePOWER")
        arrset(ielftype,ie,iet_["ePOWER"])
        vname = "ET1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="PWR",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-0.229))
        posep = findfirst(x->x=="SC",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1000.0))
        ename = "XPROD3"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD1")
        arrset(ielftype,ie,iet_["ePROD1"])
        vname = "AW2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="VA",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "IW2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="VB",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "MF2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="VC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TT2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="VD",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PW2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="VE",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "X2PF"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePOWER")
        arrset(ielftype,ie,iet_["ePOWER"])
        vname = "ET2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="PWR",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-0.146))
        posep = findfirst(x->x=="SC",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1000.0))
        ename = "X2PG"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePOWER")
        arrset(ielftype,ie,iet_["ePOWER"])
        vname = "ET2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="PWR",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.648))
        posep = findfirst(x->x=="SC",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1000.0))
        ename = "XPROD4"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD2")
        arrset(ielftype,ie,iet_["ePROD2"])
        vname = "AW2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="VA",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "MF2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="VB",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PW2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="VC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "S2L"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="VD",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "X2PL"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePOWER")
        arrset(ielftype,ie,iet_["ePOWER"])
        vname = "ET2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="PWR",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.736))
        posep = findfirst(x->x=="SC",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1000.0))
        ename = "X2PM"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePOWER")
        arrset(ielftype,ie,iet_["ePOWER"])
        vname = "ET2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="PWR",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-0.229))
        posep = findfirst(x->x=="SC",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1000.0))
        ename = "XPROD5"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD1")
        arrset(ielftype,ie,iet_["ePROD1"])
        vname = "AW3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="VA",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "IW3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="VB",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "MF3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="VC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TT3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="VD",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PW3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="VE",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "XQF"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePOWER")
        arrset(ielftype,ie,iet_["ePOWER"])
        vname = "ET3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="PWR",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.539))
        posep = findfirst(x->x=="SC",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1000.0))
        ename = "XQG"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePOWER")
        arrset(ielftype,ie,iet_["ePOWER"])
        vname = "ET3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="PWR",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.772))
        posep = findfirst(x->x=="SC",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1000.0))
        ename = "XPROD6"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD2")
        arrset(ielftype,ie,iet_["ePROD2"])
        vname = "AW3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="VA",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "MF3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="VB",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PW3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="VC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "S3L"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="VD",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "XQL"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePOWER")
        arrset(ielftype,ie,iet_["ePOWER"])
        vname = "ET3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="PWR",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-1.33))
        posep = findfirst(x->x=="SC",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(100.0))
        ename = "XQM"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePOWER")
        arrset(ielftype,ie,iet_["ePOWER"])
        vname = "ET3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="PWR",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.498))
        posep = findfirst(x->x=="SC",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(100.0))
        ename = "SMFE1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eX7Y")
        arrset(ielftype,ie,iet_["eX7Y"])
        vname = "MF1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "IW1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="Y1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PW1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="Y2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "IW2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="Y3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PW2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="Y4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "IW3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="Y5",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PW3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="Y6",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "INW"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="Y7",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "SMFE2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eX5Y")
        arrset(ielftype,ie,iet_["eX5Y"])
        vname = "MF2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "IW2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="Y1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PW2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="Y2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "IW3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="Y3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PW3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="Y4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "INW"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="Y5",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "SMFE3"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eX3Y")
        arrset(ielftype,ie,iet_["eX3Y"])
        vname = "MF3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "IW3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="Y1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PW3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="Y2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "INW"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="Y3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "TT1BT1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en2PR")
        arrset(ielftype,ie,iet_["en2PR"])
        vname = "TT1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "BT1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "TT2BT2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en2PR")
        arrset(ielftype,ie,iet_["en2PR"])
        vname = "TT2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "BT2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "TT3BT3"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en2PR")
        arrset(ielftype,ie,iet_["en2PR"])
        vname = "TT3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "BT3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "XBIG11"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eBIG1")
        arrset(ielftype,ie,iet_["eBIG1"])
        vname = "TT1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="LH",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "BT1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="TH",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PW1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="LL",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "MF1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "XBIG12"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eBIG1")
        arrset(ielftype,ie,iet_["eBIG1"])
        vname = "TT2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="LH",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "BT2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="TH",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PW2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="LL",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "MF2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "XBIG13"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eBIG1")
        arrset(ielftype,ie,iet_["eBIG1"])
        vname = "TT3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="LH",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "BT3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="TH",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PW3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="LL",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "MF3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),Float64(1.0e+4),nothing))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gSUMM",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["STA1"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["XPROD1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(5272.77))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["XPF"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(160.909))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["XPG"])
        loaset(pbm.grelw,ig,posel,Float64(282.874))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["XPROD2"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(0.64570846))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["XPL"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(31.136196))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["XPM"])
        loaset(pbm.grelw,ig,posel,Float64(12.092112))
        ig = ig_["STA2"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["XPROD3"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(5272.77))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["X2PF"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(160.909))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["X2PG"])
        loaset(pbm.grelw,ig,posel,Float64(282.874))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["XPROD4"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(0.64570846))
        ig = ig_["STA1"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["X2PL"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(31.136196))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["X2PM"])
        loaset(pbm.grelw,ig,posel,Float64(12.092112))
        ig = ig_["STA3"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["XPROD5"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(5272.77))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["XQF"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(181.806))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["XQG"])
        loaset(pbm.grelw,ig,posel,Float64(232.57))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["XPROD6"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(0.49783215))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["XQL"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-0.22424514))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["XQM"])
        loaset(pbm.grelw,ig,posel,Float64(20.708238))
        ig = ig_["LAUN"]
        arrset(pbm.grftype,ig,"gSUMM")
        ig = ig_["SMF1"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["SMFE1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["SMF2"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["SMFE2"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["SMF3"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["SMFE3"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["SI1A"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["TT1BT1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["SI1B"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["TT1BT1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["SI2A"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["TT2BT2"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["SI2B"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["TT2BT2"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["SI3A"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["TT3BT3"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["SI3B"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["TT3BT3"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["GLGCON"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["XBIG11"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-32.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["XBIG12"])
        loaset(pbm.grelw,ig,posel,Float64(-32.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["XBIG13"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-32.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[1:pb.nle] = grange[legrps]
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.clower[pb.nle+pb.neq+1:pb.m] = zeros(Float64,pb.nge)
        pb.cupper[1:pb.nge] = grange[gegrps]
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-OOR2-MY-25-28"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eBIG1"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        LG = log(EV_[4])
        f_   = (EV_[1]*EV_[2]*LG)/EV_[3]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[4] = (EV_[1]*EV_[2])/(EV_[4]*EV_[3])
            g_[1] = (EV_[2]*LG)/EV_[3]
            g_[2] = (EV_[1]*LG)/EV_[3]
            g_[3] = -(EV_[1]*EV_[2]*LG)/(EV_[3]^2)
            if nargout>2
                H_ = zeros(Float64,4,4)
                H_[4,4] = -(EV_[1]*EV_[2])/(EV_[3]*EV_[4]^2)
                H_[4,1] = EV_[2]/(EV_[4]*EV_[3])
                H_[1,4] = H_[4,1]
                H_[4,2] = EV_[1]/(EV_[4]*EV_[3])
                H_[2,4] = H_[4,2]
                H_[4,3] = -(EV_[1]*EV_[2])/(EV_[3]^2*EV_[4])
                H_[3,4] = H_[4,3]
                H_[1,2] = LG/EV_[3]
                H_[2,1] = H_[1,2]
                H_[1,3] = -EV_[2]*LG/EV_[3]^2
                H_[3,1] = H_[1,3]
                H_[2,3] = -EV_[1]*LG/EV_[3]^2
                H_[3,2] = H_[2,3]
                H_[3,3] = 2.0*(EV_[1]*EV_[2]*LG)/(EV_[3]^3.0)
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "ePROD1"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        EA = 1.2781
        VA0 = EV_[1]^EA
        VA1 = EA*EV_[1]^(EA-1.0)
        VA2 = EA*(EA-1.0)*EV_[1]^(EA-2.0)
        EB = -0.1959
        VB0 = EV_[2]^EB
        VB1 = EB*EV_[2]^(EB-1.0)
        VB2 = EB*(EB-1.0)*EV_[2]^(EB-2.0)
        EC = 2.4242
        VC0 = EV_[3]^EC
        VC1 = EC*EV_[3]^(EC-1.0)
        VC2 = EC*(EC-1.0)*EV_[3]^(EC-2.0)
        ED = 0.38745
        VD0 = EV_[4]^ED
        VD1 = ED*EV_[4]^(ED-1.0)
        VD2 = ED*(ED-1.0)*EV_[4]^(ED-2.0)
        EE = 0.9904
        VE0 = EV_[5]^EE
        VE1 = EE*EV_[5]^(EE-1.0)
        VE2 = EE*(EE-1.0)*EV_[5]^(EE-2.0)
        f_   = VA0*VB0*VC0*VD0*VE0
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = VA1*VB0*VC0*VD0*VE0
            g_[2] = VA0*VB1*VC0*VD0*VE0
            g_[3] = VA0*VB0*VC1*VD0*VE0
            g_[4] = VA0*VB0*VC0*VD1*VE0
            g_[5] = VA0*VB0*VC0*VD0*VE1
            if nargout>2
                H_ = zeros(Float64,5,5)
                H_[1,1] = VA2*VB0*VC0*VD0*VE0
                H_[1,2] = VA1*VB1*VC0*VD0*VE0
                H_[2,1] = H_[1,2]
                H_[1,3] = VA1*VB0*VC1*VD0*VE0
                H_[3,1] = H_[1,3]
                H_[1,4] = VA1*VB0*VC0*VD1*VE0
                H_[4,1] = H_[1,4]
                H_[1,5] = VA1*VB0*VC0*VD0*VE1
                H_[5,1] = H_[1,5]
                H_[2,2] = VA0*VB2*VC0*VD0*VE0
                H_[2,3] = VA0*VB1*VC1*VD0*VE0
                H_[3,2] = H_[2,3]
                H_[2,4] = VA0*VB1*VC0*VD1*VE0
                H_[4,2] = H_[2,4]
                H_[2,5] = VA0*VB1*VC0*VD0*VE1
                H_[5,2] = H_[2,5]
                H_[3,3] = VA0*VB0*VC2*VD0*VE0
                H_[3,4] = VA0*VB0*VC1*VD1*VE0
                H_[4,3] = H_[3,4]
                H_[3,5] = VA0*VB0*VC1*VD0*VE1
                H_[5,3] = H_[3,5]
                H_[4,4] = VA0*VB0*VC0*VD2*VE0
                H_[4,5] = VA0*VB0*VC0*VD1*VE1
                H_[5,4] = H_[4,5]
                H_[5,5] = VA0*VB0*VC0*VD0*VE2
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "ePROD2"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        EA = 0.3322
        VA0 = EV_[1]^EA
        VA1 = EA*EV_[1]^(EA-1.0)
        VA2 = EA*(EA-1.0)*EV_[1]^(EA-2.0)
        EB = -1.5935
        VB0 = EV_[2]^EB
        VB1 = EB*EV_[2]^(EB-1.0)
        VB2 = EB*(EB-1.0)*EV_[2]^(EB-2.0)
        EC = 0.2363
        VC0 = EV_[3]^EC
        VC1 = EC*EV_[3]^(EC-1.0)
        VC2 = EC*(EC-1.0)*EV_[3]^(EC-2.0)
        ED = 0.1079
        VD0 = EV_[4]^ED
        VD1 = ED*EV_[4]^(ED-1.0)
        VD2 = ED*(ED-1.0)*EV_[4]^(ED-2.0)
        f_   = VA0*VB0*VC0*VD0
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = VA1*VB0*VC0*VD0
            g_[2] = VA0*VB1*VC0*VD0
            g_[3] = VA0*VB0*VC1*VD0
            g_[4] = VA0*VB0*VC0*VD1
            if nargout>2
                H_ = zeros(Float64,4,4)
                H_[1,1] = VA2*VB0*VC0*VD0
                H_[1,2] = VA1*VB1*VC0*VD0
                H_[2,1] = H_[1,2]
                H_[1,3] = VA1*VB0*VC1*VD0
                H_[3,1] = H_[1,3]
                H_[1,4] = VA1*VB0*VC0*VD1
                H_[4,1] = H_[1,4]
                H_[2,2] = VA0*VB2*VC0*VD0
                H_[2,3] = VA0*VB1*VC1*VD0
                H_[3,2] = H_[2,3]
                H_[2,4] = VA0*VB1*VC0*VD1
                H_[4,2] = H_[2,4]
                H_[3,3] = VA0*VB0*VC2*VD0
                H_[3,4] = VA0*VB0*VC1*VD1
                H_[4,3] = H_[3,4]
                H_[4,4] = VA0*VB0*VC0*VD2
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "ePOWER"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        SCPWR = pbm.elpar[iel_][1]/(pbm.elpar[iel_][2]^pbm.elpar[iel_][1])
        f_   = (EV_[1]/pbm.elpar[iel_][2])^pbm.elpar[iel_][1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = SCPWR*EV_[1]^(pbm.elpar[iel_][1]-1.0)
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = SCPWR*(pbm.elpar[iel_][1]-1.0)*EV_[1]^(pbm.elpar[iel_][1]-2.0)
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eX7Y"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,2,8)
        IV_ =  zeros(Float64,2)
        U_[1,1] = U_[1,1]+1
        U_[2,2] = U_[2,2]+1
        U_[2,3] = U_[2,3]+1
        U_[2,4] = U_[2,4]+1
        U_[2,5] = U_[2,5]+1
        U_[2,6] = U_[2,6]+1
        U_[2,7] = U_[2,7]+1
        U_[2,8] = U_[2,8]+1
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        f_   = IV_[1]*IV_[2]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = IV_[2]
            g_[2] = IV_[1]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = 1.0
                H_[2,1] = H_[1,2]
                H_ = U_'*H_*U_
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eX5Y"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,2,6)
        IV_ =  zeros(Float64,2)
        U_[1,1] = U_[1,1]+1
        U_[2,2] = U_[2,2]+1
        U_[2,3] = U_[2,3]+1
        U_[2,4] = U_[2,4]+1
        U_[2,5] = U_[2,5]+1
        U_[2,6] = U_[2,6]+1
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        f_   = IV_[1]*IV_[2]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = IV_[2]
            g_[2] = IV_[1]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = 1.0
                H_[2,1] = H_[1,2]
                H_ = U_'*H_*U_
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eX3Y"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,2,4)
        IV_ =  zeros(Float64,2)
        U_[1,1] = U_[1,1]+1
        U_[2,2] = U_[2,2]+1
        U_[2,3] = U_[2,3]+1
        U_[2,4] = U_[2,4]+1
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        f_   = IV_[1]*IV_[2]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = IV_[2]
            g_[2] = IV_[1]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = 1.0
                H_[2,1] = H_[1,2]
                H_ = U_'*H_*U_
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "en2PR"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[2]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[2]
            g_[2] = EV_[1]
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = 1.0
                H_[2,1] = H_[1,2]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    elseif action == "gSUMM"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_= GVAR_^0.460
        if nargout>1
            g_ = 0.460*GVAR_^(-0.540)
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = -0.2484*GVAR_^(-1.540)
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

