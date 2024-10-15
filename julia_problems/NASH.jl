function NASH(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : NASH
#    *********
# 
#    A quadratic programming reformulation of a linear
#    complementarity problem arising from Nash equilibrium
#    provided by Michael Ferris
# 
#    classification = "C-QLR2-AN-72-24"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "NASH"

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
        v_["N"] = 72
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for J = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("X"*string(J),ix_)
            arrset(pb.xnames,iv,"X"*string(J))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X56"]
        pbm.A[ig,iv] += Float64(1000.0)
        iv = ix_["X57"]
        pbm.A[ig,iv] += Float64(500.0)
        iv = ix_["X58"]
        pbm.A[ig,iv] += Float64(1000.0)
        ig,ig_,_ = s2mpj_ii("C1",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C1")
        iv = ix_["X25"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X49"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X3"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X5"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X6"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X7"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("C2",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C2")
        iv = ix_["X26"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X50"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(0.02309)
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(0.02309)
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(0.02309)
        iv = ix_["X6"]
        pbm.A[ig,iv] += Float64(0.02309)
        iv = ix_["X11"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X22"]
        pbm.A[ig,iv] += Float64(0.288626)
        iv = ix_["X23"]
        pbm.A[ig,iv] += Float64(0.263887)
        iv = ix_["X24"]
        pbm.A[ig,iv] += Float64(0.447486)
        ig,ig_,_ = s2mpj_ii("C3",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C3")
        iv = ix_["X27"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X51"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(0.02309)
        iv = ix_["X3"]
        pbm.A[ig,iv] += Float64(0.02309)
        iv = ix_["X5"]
        pbm.A[ig,iv] += Float64(0.02309)
        iv = ix_["X7"]
        pbm.A[ig,iv] += Float64(0.02309)
        iv = ix_["X12"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X22"]
        pbm.A[ig,iv] += Float64(0.288626)
        iv = ix_["X23"]
        pbm.A[ig,iv] += Float64(0.263887)
        iv = ix_["X24"]
        pbm.A[ig,iv] += Float64(0.447486)
        ig,ig_,_ = s2mpj_ii("C4",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C4")
        iv = ix_["X28"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X52"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(0.02309)
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(0.02309)
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(0.02309)
        iv = ix_["X6"]
        pbm.A[ig,iv] += Float64(0.02309)
        iv = ix_["X11"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X22"]
        pbm.A[ig,iv] += Float64(0.288626)
        iv = ix_["X23"]
        pbm.A[ig,iv] += Float64(0.263887)
        iv = ix_["X24"]
        pbm.A[ig,iv] += Float64(0.447486)
        ig,ig_,_ = s2mpj_ii("C5",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C5")
        iv = ix_["X29"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X53"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(0.02309)
        iv = ix_["X3"]
        pbm.A[ig,iv] += Float64(0.02309)
        iv = ix_["X5"]
        pbm.A[ig,iv] += Float64(0.02309)
        iv = ix_["X7"]
        pbm.A[ig,iv] += Float64(0.02309)
        iv = ix_["X12"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X22"]
        pbm.A[ig,iv] += Float64(0.288626)
        iv = ix_["X23"]
        pbm.A[ig,iv] += Float64(0.263887)
        iv = ix_["X24"]
        pbm.A[ig,iv] += Float64(0.447486)
        ig,ig_,_ = s2mpj_ii("C6",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C6")
        iv = ix_["X30"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X54"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(0.02309)
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(0.02309)
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(0.02309)
        iv = ix_["X6"]
        pbm.A[ig,iv] += Float64(0.02309)
        iv = ix_["X11"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X22"]
        pbm.A[ig,iv] += Float64(0.288626)
        iv = ix_["X23"]
        pbm.A[ig,iv] += Float64(0.263887)
        iv = ix_["X24"]
        pbm.A[ig,iv] += Float64(0.447486)
        ig,ig_,_ = s2mpj_ii("C7",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C7")
        iv = ix_["X31"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X55"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(0.02309)
        iv = ix_["X3"]
        pbm.A[ig,iv] += Float64(0.02309)
        iv = ix_["X5"]
        pbm.A[ig,iv] += Float64(0.02309)
        iv = ix_["X7"]
        pbm.A[ig,iv] += Float64(0.02309)
        iv = ix_["X12"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X22"]
        pbm.A[ig,iv] += Float64(0.288626)
        iv = ix_["X23"]
        pbm.A[ig,iv] += Float64(0.263887)
        iv = ix_["X24"]
        pbm.A[ig,iv] += Float64(0.447486)
        ig,ig_,_ = s2mpj_ii("C8",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C8")
        iv = ix_["X32"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X56"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X11"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X22"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("C9",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C9")
        iv = ix_["X33"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X57"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X11"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X23"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("C10",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C10")
        iv = ix_["X34"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X58"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X12"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X23"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("C11",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C11")
        iv = ix_["X35"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X59"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X6"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X8"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X9"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("C12",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C12")
        iv = ix_["X36"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X60"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X3"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X5"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X7"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X10"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("C13",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C13")
        iv = ix_["X37"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X61"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X16"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X19"]
        pbm.A[ig,iv] += Float64(-.33)
        iv = ix_["X20"]
        pbm.A[ig,iv] += Float64(0.33)
        ig,ig_,_ = s2mpj_ii("C14",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C14")
        iv = ix_["X38"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X62"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X17"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X19"]
        pbm.A[ig,iv] += Float64(-.67)
        iv = ix_["X20"]
        pbm.A[ig,iv] += Float64(-.33)
        ig,ig_,_ = s2mpj_ii("C15",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C15")
        iv = ix_["X39"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X63"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X18"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X19"]
        pbm.A[ig,iv] += Float64(-.33)
        iv = ix_["X20"]
        pbm.A[ig,iv] += Float64(-.67)
        ig,ig_,_ = s2mpj_ii("C16",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C16")
        iv = ix_["X40"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X64"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X13"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("C17",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C17")
        iv = ix_["X41"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X65"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X14"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("C18",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C18")
        iv = ix_["X42"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X66"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X15"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("C19",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C19")
        iv = ix_["X43"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X67"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X13"]
        pbm.A[ig,iv] += Float64(0.33)
        iv = ix_["X14"]
        pbm.A[ig,iv] += Float64(0.67)
        iv = ix_["X15"]
        pbm.A[ig,iv] += Float64(0.33)
        iv = ix_["X22"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("C20",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C20")
        iv = ix_["X44"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X68"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X13"]
        pbm.A[ig,iv] += Float64(-.33)
        iv = ix_["X14"]
        pbm.A[ig,iv] += Float64(0.33)
        iv = ix_["X15"]
        pbm.A[ig,iv] += Float64(0.67)
        iv = ix_["X23"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("C21",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C21")
        iv = ix_["X45"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X69"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X24"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("C22",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C22")
        iv = ix_["X46"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X70"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(-.288626)
        iv = ix_["X8"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X19"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X22"]
        pbm.A[ig,iv] += Float64(8.892169)
        iv = ix_["X23"]
        pbm.A[ig,iv] += Float64(-3.298588)
        iv = ix_["X24"]
        pbm.A[ig,iv] += Float64(-5.593581)
        ig,ig_,_ = s2mpj_ii("C23",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C23")
        iv = ix_["X47"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X71"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(-.263887)
        iv = ix_["X9"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X10"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X20"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X22"]
        pbm.A[ig,iv] += Float64(-3.298588)
        iv = ix_["X23"]
        pbm.A[ig,iv] += Float64(8.412719)
        iv = ix_["X24"]
        pbm.A[ig,iv] += Float64(-5.114131)
        ig,ig_,_ = s2mpj_ii("C24",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C24")
        iv = ix_["X48"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X72"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(-.447486)
        iv = ix_["X21"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X22"]
        pbm.A[ig,iv] += Float64(-5.593581)
        iv = ix_["X23"]
        pbm.A[ig,iv] += Float64(-5.114131)
        iv = ix_["X24"]
        pbm.A[ig,iv] += Float64(10.707712)
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
        pbm.gconst[ig_["C2"]] = Float64(35.100673)
        pbm.gconst[ig_["C3"]] = Float64(35.100673)
        pbm.gconst[ig_["C4"]] = Float64(35.100673)
        pbm.gconst[ig_["C5"]] = Float64(35.100673)
        pbm.gconst[ig_["C6"]] = Float64(35.100673)
        pbm.gconst[ig_["C7"]] = Float64(35.100673)
        pbm.gconst[ig_["C8"]] = Float64(-15.0)
        pbm.gconst[ig_["C9"]] = Float64(-15.0)
        pbm.gconst[ig_["C10"]] = Float64(-20.0)
        pbm.gconst[ig_["C22"]] = Float64(61.241589)
        pbm.gconst[ig_["C23"]] = Float64(-1.150548)
        pbm.gconst[ig_["C24"]] = Float64(-60.091041)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = fill(0.E+00,pb.n)
        pb.xupper = fill(0.E+00,pb.n)
        pb.xlower[ix_["X1"]] = -Inf
        pb.xupper[ix_["X1"]] = +Inf
        pb.xupper[ix_["X8"]] = 1000.0
        pb.xupper[ix_["X9"]] = 500.0
        pb.xupper[ix_["X10"]] = 1000.0
        pb.xlower[ix_["X11"]] = -Inf
        pb.xupper[ix_["X11"]] = +Inf
        pb.xlower[ix_["X12"]] = -Inf
        pb.xupper[ix_["X12"]] = +Inf
        pb.xlower[ix_["X13"]] = -Inf
        pb.xupper[ix_["X13"]] = +Inf
        pb.xlower[ix_["X14"]] = -Inf
        pb.xupper[ix_["X14"]] = +Inf
        pb.xlower[ix_["X15"]] = -Inf
        pb.xupper[ix_["X15"]] = +Inf
        pb.xlower[ix_["X16"]] = -Inf
        pb.xupper[ix_["X16"]] = +Inf
        pb.xlower[ix_["X17"]] = -Inf
        pb.xupper[ix_["X17"]] = +Inf
        pb.xlower[ix_["X18"]] = -Inf
        pb.xupper[ix_["X18"]] = +Inf
        pb.xlower[ix_["X19"]] = -Inf
        pb.xupper[ix_["X19"]] = +Inf
        pb.xlower[ix_["X20"]] = -Inf
        pb.xupper[ix_["X20"]] = +Inf
        pb.xlower[ix_["X21"]] = -Inf
        pb.xupper[ix_["X21"]] = +Inf
        pb.xlower[ix_["X22"]] = -Inf
        pb.xupper[ix_["X22"]] = +Inf
        pb.xlower[ix_["X23"]] = -Inf
        pb.xupper[ix_["X23"]] = +Inf
        pb.xlower[ix_["X24"]] = -Inf
        pb.xupper[ix_["X24"]] = +Inf
        #%%%%%%%%%%%%%%%%%%%% QUADRATIC %%%%%%%%%%%%%%%%%%%
        ix1 = ix_["X25"]
        ix2 = ix_["X1"]
        pbm.H[ix1,ix2] = Float64(1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X26"]
        ix2 = ix_["X2"]
        pbm.H[ix1,ix2] = Float64(1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X27"]
        ix2 = ix_["X3"]
        pbm.H[ix1,ix2] = Float64(1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X28"]
        ix2 = ix_["X4"]
        pbm.H[ix1,ix2] = Float64(1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X29"]
        ix2 = ix_["X5"]
        pbm.H[ix1,ix2] = Float64(1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X30"]
        ix2 = ix_["X6"]
        pbm.H[ix1,ix2] = Float64(1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X31"]
        ix2 = ix_["X7"]
        pbm.H[ix1,ix2] = Float64(1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X32"]
        ix2 = ix_["X8"]
        pbm.H[ix1,ix2] = Float64(1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X33"]
        ix2 = ix_["X9"]
        pbm.H[ix1,ix2] = Float64(1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X34"]
        ix2 = ix_["X10"]
        pbm.H[ix1,ix2] = Float64(1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X35"]
        ix2 = ix_["X11"]
        pbm.H[ix1,ix2] = Float64(1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X36"]
        ix2 = ix_["X12"]
        pbm.H[ix1,ix2] = Float64(1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X37"]
        ix2 = ix_["X13"]
        pbm.H[ix1,ix2] = Float64(1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X38"]
        ix2 = ix_["X14"]
        pbm.H[ix1,ix2] = Float64(1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X39"]
        ix2 = ix_["X15"]
        pbm.H[ix1,ix2] = Float64(1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X40"]
        ix2 = ix_["X16"]
        pbm.H[ix1,ix2] = Float64(1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X41"]
        ix2 = ix_["X17"]
        pbm.H[ix1,ix2] = Float64(1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X42"]
        ix2 = ix_["X18"]
        pbm.H[ix1,ix2] = Float64(1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X43"]
        ix2 = ix_["X19"]
        pbm.H[ix1,ix2] = Float64(1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X44"]
        ix2 = ix_["X20"]
        pbm.H[ix1,ix2] = Float64(1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X45"]
        ix2 = ix_["X21"]
        pbm.H[ix1,ix2] = Float64(1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X46"]
        ix2 = ix_["X22"]
        pbm.H[ix1,ix2] = Float64(1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X47"]
        ix2 = ix_["X23"]
        pbm.H[ix1,ix2] = Float64(1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X48"]
        ix2 = ix_["X24"]
        pbm.H[ix1,ix2] = Float64(1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X49"]
        ix2 = ix_["X1"]
        pbm.H[ix1,ix2] = Float64(-1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X50"]
        ix2 = ix_["X2"]
        pbm.H[ix1,ix2] = Float64(-1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X51"]
        ix2 = ix_["X3"]
        pbm.H[ix1,ix2] = Float64(-1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X52"]
        ix2 = ix_["X4"]
        pbm.H[ix1,ix2] = Float64(-1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X53"]
        ix2 = ix_["X5"]
        pbm.H[ix1,ix2] = Float64(-1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X54"]
        ix2 = ix_["X6"]
        pbm.H[ix1,ix2] = Float64(-1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X55"]
        ix2 = ix_["X7"]
        pbm.H[ix1,ix2] = Float64(-1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X56"]
        ix2 = ix_["X8"]
        pbm.H[ix1,ix2] = Float64(-1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X57"]
        ix2 = ix_["X9"]
        pbm.H[ix1,ix2] = Float64(-1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X58"]
        ix2 = ix_["X10"]
        pbm.H[ix1,ix2] = Float64(-1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X59"]
        ix2 = ix_["X11"]
        pbm.H[ix1,ix2] = Float64(-1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X60"]
        ix2 = ix_["X12"]
        pbm.H[ix1,ix2] = Float64(-1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X61"]
        ix2 = ix_["X13"]
        pbm.H[ix1,ix2] = Float64(-1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X62"]
        ix2 = ix_["X14"]
        pbm.H[ix1,ix2] = Float64(-1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X63"]
        ix2 = ix_["X15"]
        pbm.H[ix1,ix2] = Float64(-1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X64"]
        ix2 = ix_["X16"]
        pbm.H[ix1,ix2] = Float64(-1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X65"]
        ix2 = ix_["X17"]
        pbm.H[ix1,ix2] = Float64(-1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X66"]
        ix2 = ix_["X18"]
        pbm.H[ix1,ix2] = Float64(-1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X67"]
        ix2 = ix_["X19"]
        pbm.H[ix1,ix2] = Float64(-1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X68"]
        ix2 = ix_["X20"]
        pbm.H[ix1,ix2] = Float64(-1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X69"]
        ix2 = ix_["X21"]
        pbm.H[ix1,ix2] = Float64(-1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X70"]
        ix2 = ix_["X22"]
        pbm.H[ix1,ix2] = Float64(-1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X71"]
        ix2 = ix_["X23"]
        pbm.H[ix1,ix2] = Float64(-1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_["X72"]
        ix2 = ix_["X24"]
        pbm.H[ix1,ix2] = Float64(-1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        Hsave = pbm.H[ 1:pb.n, 1:pb.n ]
        pbm.H = Hsave
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = collect(1:length(pbm.congrps))
        pb.pbclass = "C-QLR2-AN-72-24"
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

