function DALLASS(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DALLASS
#    *********
# 
#    The small Dallas water distribution problem
#    The problem is also named "W30" in some references.
#    This is a nonlinear network problem with conditioning of
#    the order of 10**4.
# 
#    Source:
#    R. Dembo,
#    private communication, 1986.
# 
#    SIF input: Ph. Toint, June 1990.
# 
#    classification = "C-ONR2-MN-46-31"
# 
#    Number of arcs
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "DALLASS"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["N"] = 46
        v_["NODES"] = 31
        v_["1"] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X42"]
        pbm.A[ig,iv] += Float64(-6.38400e+02)
        iv = ix_["X43"]
        pbm.A[ig,iv] += Float64(-6.33000e+02)
        iv = ix_["X44"]
        pbm.A[ig,iv] += Float64(-5.54500e+02)
        iv = ix_["X45"]
        pbm.A[ig,iv] += Float64(-5.05000e+02)
        iv = ix_["X46"]
        pbm.A[ig,iv] += Float64(-4.36900e+02)
        ig,ig_,_ = s2mpj_ii("N1",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"N1")
        iv = ix_["X46"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X41"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("N2",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"N2")
        iv = ix_["X45"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("N3",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"N3")
        iv = ix_["X44"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X3"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("N4",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"N4")
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("N5",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"N5")
        iv = ix_["X16"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X7"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X6"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X5"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("N6",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"N6")
        iv = ix_["X7"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X9"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X8"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("N7",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"N7")
        iv = ix_["X9"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X10"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("N8",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"N8")
        iv = ix_["X10"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X12"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X11"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("N9",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"N9")
        iv = ix_["X12"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X3"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X13"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("N10",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"N10")
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X16"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X15"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X14"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("N11",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"N11")
        iv = ix_["X15"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X13"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X5"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X17"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("N12",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"N12")
        iv = ix_["X20"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X19"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X8"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X6"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X18"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("N13",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"N13")
        iv = ix_["X42"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X18"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X19"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("N14",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"N14")
        iv = ix_["X21"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X20"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("N15",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"N15")
        iv = ix_["X43"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X21"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("N16",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"N16")
        iv = ix_["X14"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X11"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X23"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X22"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("N17",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"N17")
        iv = ix_["X23"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X25"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X24"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("N18",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"N18")
        iv = ix_["X31"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X25"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X22"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X26"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("N19",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"N19")
        iv = ix_["X26"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X17"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X28"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X27"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("N20",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"N20")
        iv = ix_["X28"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("N21",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"N21")
        iv = ix_["X31"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X30"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X29"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("N22",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"N22")
        iv = ix_["X30"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X27"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("N23",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"N23")
        iv = ix_["X24"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X32"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("N24",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"N24")
        iv = ix_["X38"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X29"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X34"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X33"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("N25",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"N25")
        iv = ix_["X32"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X35"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("N26",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"N26")
        iv = ix_["X35"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X37"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X36"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("N27",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"N27")
        iv = ix_["X37"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X34"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("N28",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"N28")
        iv = ix_["X36"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X40"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X39"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X38"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("N29",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"N29")
        iv = ix_["X39"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X33"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("N30",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"N30")
        iv = ix_["X40"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X41"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("N31",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"N31")
        iv = ix_["X46"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X45"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X44"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X43"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X42"]
        pbm.A[ig,iv] += Float64(-1.0)
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
        pbm.gconst[ig_["N5"]] = Float64(2.80000)
        pbm.gconst[ig_["N7"]] = Float64(4.03000e-01)
        pbm.gconst[ig_["N8"]] = Float64(5.92000e-01)
        pbm.gconst[ig_["N9"]] = Float64(1.15600)
        pbm.gconst[ig_["N10"]] = Float64(2.00000e-01)
        pbm.gconst[ig_["N11"]] = Float64(4.95000e-01)
        pbm.gconst[ig_["N16"]] = Float64(3.13000e-01)
        pbm.gconst[ig_["N17"]] = Float64(8.44000e-01)
        pbm.gconst[ig_["N18"]] = Float64(3.31000e-01)
        pbm.gconst[ig_["N19"]] = Float64(5.30000e-02)
        pbm.gconst[ig_["N21"]] = Float64(2.72000e-01)
        pbm.gconst[ig_["N22"]] = Float64(8.83000e-01)
        pbm.gconst[ig_["N23"]] = Float64(5.71000e-01)
        pbm.gconst[ig_["N24"]] = Float64(7.55000e-01)
        pbm.gconst[ig_["N26"]] = Float64(5.27000e-01)
        pbm.gconst[ig_["N29"]] = Float64(1.00000e-03)
        pbm.gconst[ig_["N31"]] = Float64(-1.01960e+01)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = fill(-2.00000e+02,pb.n)
        pb.xupper = fill(2.00000e+02,pb.n)
        pb.xlower[ix_["X1"]] = 0.00000
        pb.xupper[ix_["X1"]] = 2.11673e+01
        pb.xlower[ix_["X2"]] = 0.00000
        pb.xupper[ix_["X2"]] = 4.37635e+01
        pb.xlower[ix_["X3"]] = 0.00000
        pb.xupper[ix_["X3"]] = 3.28255e+01
        pb.xlower[ix_["X19"]] = 0.00000
        pb.xupper[ix_["X19"]] = 2.20120e+01
        pb.xlower[ix_["X21"]] = 0.00000
        pb.xupper[ix_["X21"]] = 1.36703e+01
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(-2.00000e+02),pb.n)
        if haskey(ix_,"X1")
            pb.x0[ix_["X1"]] = Float64(2.11673e+01)
        else
            pb.y0[findfirst(x->x==ig_["X1"],pbm.congrps)] = Float64(2.11673e+01)
        end
        if haskey(ix_,"X2")
            pb.x0[ix_["X2"]] = Float64(4.37635e+01)
        else
            pb.y0[findfirst(x->x==ig_["X2"],pbm.congrps)] = Float64(4.37635e+01)
        end
        if haskey(ix_,"X3")
            pb.x0[ix_["X3"]] = Float64(3.28255e+01)
        else
            pb.y0[findfirst(x->x==ig_["X3"],pbm.congrps)] = Float64(3.28255e+01)
        end
        if haskey(ix_,"X4")
            pb.x0[ix_["X4"]] = Float64(1.42109e-14)
        else
            pb.y0[findfirst(x->x==ig_["X4"],pbm.congrps)] = Float64(1.42109e-14)
        end
        if haskey(ix_,"X5")
            pb.x0[ix_["X5"]] = Float64(1.68826e+02)
        else
            pb.y0[findfirst(x->x==ig_["X5"],pbm.congrps)] = Float64(1.68826e+02)
        end
        if haskey(ix_,"X7")
            pb.x0[ix_["X7"]] = Float64(2.81745e+01)
        else
            pb.y0[findfirst(x->x==ig_["X7"],pbm.congrps)] = Float64(2.81745e+01)
        end
        if haskey(ix_,"X8")
            pb.x0[ix_["X8"]] = Float64(8.75603e+01)
        else
            pb.y0[findfirst(x->x==ig_["X8"],pbm.congrps)] = Float64(8.75603e+01)
        end
        if haskey(ix_,"X9")
            pb.x0[ix_["X9"]] = Float64(-5.93858e+01)
        else
            pb.y0[findfirst(x->x==ig_["X9"],pbm.congrps)] = Float64(-5.93858e+01)
        end
        if haskey(ix_,"X10")
            pb.x0[ix_["X10"]] = Float64(-5.97888e+01)
        else
            pb.y0[findfirst(x->x==ig_["X10"],pbm.congrps)] = Float64(-5.97888e+01)
        end
        if haskey(ix_,"X11")
            pb.x0[ix_["X11"]] = Float64(1.83383e+02)
        else
            pb.y0[findfirst(x->x==ig_["X11"],pbm.congrps)] = Float64(1.83383e+02)
        end
        if haskey(ix_,"X13")
            pb.x0[ix_["X13"]] = Float64(-1.68331e+02)
        else
            pb.y0[findfirst(x->x==ig_["X13"],pbm.congrps)] = Float64(-1.68331e+02)
        end
        if haskey(ix_,"X15")
            pb.x0[ix_["X15"]] = Float64(2.00000e+02)
        else
            pb.y0[findfirst(x->x==ig_["X15"],pbm.congrps)] = Float64(2.00000e+02)
        end
        if haskey(ix_,"X16")
            pb.x0[ix_["X16"]] = Float64(2.00000e-01)
        else
            pb.y0[findfirst(x->x==ig_["X16"],pbm.congrps)] = Float64(2.00000e-01)
        end
        if haskey(ix_,"X17")
            pb.x0[ix_["X17"]] = Float64(2.00000e+02)
        else
            pb.y0[findfirst(x->x==ig_["X17"],pbm.congrps)] = Float64(2.00000e+02)
        end
        if haskey(ix_,"X18")
            pb.x0[ix_["X18"]] = Float64(-7.67574e+01)
        else
            pb.y0[findfirst(x->x==ig_["X18"],pbm.congrps)] = Float64(-7.67574e+01)
        end
        if haskey(ix_,"X19")
            pb.x0[ix_["X19"]] = Float64(2.20120e+01)
        else
            pb.y0[findfirst(x->x==ig_["X19"],pbm.congrps)] = Float64(2.20120e+01)
        end
        if haskey(ix_,"X20")
            pb.x0[ix_["X20"]] = Float64(1.36703e+01)
        else
            pb.y0[findfirst(x->x==ig_["X20"],pbm.congrps)] = Float64(1.36703e+01)
        end
        if haskey(ix_,"X21")
            pb.x0[ix_["X21"]] = Float64(1.36703e+01)
        else
            pb.y0[findfirst(x->x==ig_["X21"],pbm.congrps)] = Float64(1.36703e+01)
        end
        if haskey(ix_,"X22")
            pb.x0[ix_["X22"]] = Float64(-1.98461e+02)
        else
            pb.y0[findfirst(x->x==ig_["X22"],pbm.congrps)] = Float64(-1.98461e+02)
        end
        if haskey(ix_,"X23")
            pb.x0[ix_["X23"]] = Float64(1.81531e+02)
        else
            pb.y0[findfirst(x->x==ig_["X23"],pbm.congrps)] = Float64(1.81531e+02)
        end
        if haskey(ix_,"X24")
            pb.x0[ix_["X24"]] = Float64(-1.93133e+01)
        else
            pb.y0[findfirst(x->x==ig_["X24"],pbm.congrps)] = Float64(-1.93133e+01)
        end
        if haskey(ix_,"X25")
            pb.x0[ix_["X25"]] = Float64(2.00000e+02)
        else
            pb.y0[findfirst(x->x==ig_["X25"],pbm.congrps)] = Float64(2.00000e+02)
        end
        if haskey(ix_,"X26")
            pb.x0[ix_["X26"]] = Float64(-1.98792e+02)
        else
            pb.y0[findfirst(x->x==ig_["X26"],pbm.congrps)] = Float64(-1.98792e+02)
        end
        if haskey(ix_,"X27")
            pb.x0[ix_["X27"]] = Float64(1.15500)
        else
            pb.y0[findfirst(x->x==ig_["X27"],pbm.congrps)] = Float64(1.15500)
        end
        if haskey(ix_,"X28")
            pb.x0[ix_["X28"]] = Float64(0.00000)
        else
            pb.y0[findfirst(x->x==ig_["X28"],pbm.congrps)] = Float64(0.00000)
        end
        if haskey(ix_,"X29")
            pb.x0[ix_["X29"]] = Float64(2.00000e+02)
        else
            pb.y0[findfirst(x->x==ig_["X29"],pbm.congrps)] = Float64(2.00000e+02)
        end
        if haskey(ix_,"X30")
            pb.x0[ix_["X30"]] = Float64(2.72000e-01)
        else
            pb.y0[findfirst(x->x==ig_["X30"],pbm.congrps)] = Float64(2.72000e-01)
        end
        if haskey(ix_,"X32")
            pb.x0[ix_["X32"]] = Float64(-1.98843e+01)
        else
            pb.y0[findfirst(x->x==ig_["X32"],pbm.congrps)] = Float64(-1.98843e+01)
        end
        if haskey(ix_,"X33")
            pb.x0[ix_["X33"]] = Float64(1.78834e+02)
        else
            pb.y0[findfirst(x->x==ig_["X33"],pbm.congrps)] = Float64(1.78834e+02)
        end
        if haskey(ix_,"X34")
            pb.x0[ix_["X34"]] = Float64(-1.79589e+02)
        else
            pb.y0[findfirst(x->x==ig_["X34"],pbm.congrps)] = Float64(-1.79589e+02)
        end
        if haskey(ix_,"X35")
            pb.x0[ix_["X35"]] = Float64(-1.98843e+01)
        else
            pb.y0[findfirst(x->x==ig_["X35"],pbm.congrps)] = Float64(-1.98843e+01)
        end
        if haskey(ix_,"X37")
            pb.x0[ix_["X37"]] = Float64(1.79589e+02)
        else
            pb.y0[findfirst(x->x==ig_["X37"],pbm.congrps)] = Float64(1.79589e+02)
        end
        if haskey(ix_,"X40")
            pb.x0[ix_["X40"]] = Float64(2.00000e+02)
        else
            pb.y0[findfirst(x->x==ig_["X40"],pbm.congrps)] = Float64(2.00000e+02)
        end
        if haskey(ix_,"X41")
            pb.x0[ix_["X41"]] = Float64(2.00000e+02)
        else
            pb.y0[findfirst(x->x==ig_["X41"],pbm.congrps)] = Float64(2.00000e+02)
        end
        if haskey(ix_,"X42")
            pb.x0[ix_["X42"]] = Float64(9.87694e+01)
        else
            pb.y0[findfirst(x->x==ig_["X42"],pbm.congrps)] = Float64(9.87694e+01)
        end
        if haskey(ix_,"X43")
            pb.x0[ix_["X43"]] = Float64(1.36703e+01)
        else
            pb.y0[findfirst(x->x==ig_["X43"],pbm.congrps)] = Float64(1.36703e+01)
        end
        if haskey(ix_,"X44")
            pb.x0[ix_["X44"]] = Float64(3.28255e+01)
        else
            pb.y0[findfirst(x->x==ig_["X44"],pbm.congrps)] = Float64(3.28255e+01)
        end
        if haskey(ix_,"X45")
            pb.x0[ix_["X45"]] = Float64(4.37635e+01)
        else
            pb.y0[findfirst(x->x==ig_["X45"],pbm.congrps)] = Float64(4.37635e+01)
        end
        if haskey(ix_,"X46")
            pb.x0[ix_["X46"]] = Float64(-1.78833e+02)
        else
            pb.y0[findfirst(x->x==ig_["X46"],pbm.congrps)] = Float64(-1.78833e+02)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eT1", iet_)
        loaset(elftv,it,1,"ARC")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"C1")
        loaset(elftp,it,2,"C2")
        loaset(elftp,it,3,"C3")
        it,iet_,_ = s2mpj_ii( "eT4", iet_)
        loaset(elftv,it,1,"ARC")
        loaset(elftp,it,1,"C1")
        loaset(elftp,it,2,"C2")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "E1"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eT4")
        arrset(ielftype,ie,iet_["eT4"])
        vname = "X1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(4.48060e+02))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(2.51200e+02))
        ename = "E2"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eT4")
        arrset(ielftype,ie,iet_["eT4"])
        vname = "X2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.91526e+03))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(6.46300e+01))
        ename = "E3"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eT4")
        arrset(ielftype,ie,iet_["eT4"])
        vname = "X3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.07752e+03))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(4.81400e+01))
        ename = "E4"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eT1")
            arrset(ielftype,ie,iet_["eT1"])
        end
        vname = "X4"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(6.90000e+03))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(6.60000e+01))
        posep = findfirst(x->x=="C3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.22000e+02))
        ename = "E5"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eT1")
            arrset(ielftype,ie,iet_["eT1"])
        end
        vname = "X5"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.00000e+02))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.00000e+01))
        posep = findfirst(x->x=="C3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.00000e+02))
        ename = "E6"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eT1")
            arrset(ielftype,ie,iet_["eT1"])
        end
        vname = "X6"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(6.63000e+02))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(3.60000e+01))
        posep = findfirst(x->x=="C3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.20000e+02))
        ename = "E7"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eT1")
            arrset(ielftype,ie,iet_["eT1"])
        end
        vname = "X7"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(5.10000e+03))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(6.60000e+01))
        posep = findfirst(x->x=="C3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.22000e+02))
        ename = "E8"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eT1")
            arrset(ielftype,ie,iet_["eT1"])
        end
        vname = "X8"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(6.45000e+02))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(3.00000e+01))
        posep = findfirst(x->x=="C3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.20000e+02))
        ename = "E9"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eT1")
            arrset(ielftype,ie,iet_["eT1"])
        end
        vname = "X9"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(7.40000e+03))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(6.60000e+01))
        posep = findfirst(x->x=="C3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.22000e+02))
        ename = "E10"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eT1")
            arrset(ielftype,ie,iet_["eT1"])
        end
        vname = "X10"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(5.00000e+03))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(6.60000e+01))
        posep = findfirst(x->x=="C3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(9.50000e+01))
        ename = "E11"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eT1")
            arrset(ielftype,ie,iet_["eT1"])
        end
        vname = "X11"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(8.00000e+02))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(5.40000e+01))
        posep = findfirst(x->x=="C3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.07000e+02))
        ename = "E12"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eT1")
            arrset(ielftype,ie,iet_["eT1"])
        end
        vname = "X12"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(5.20000e+03))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(4.80000e+01))
        posep = findfirst(x->x=="C3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.10000e+02))
        ename = "E13"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eT1")
            arrset(ielftype,ie,iet_["eT1"])
        end
        vname = "X13"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(6.00000e+03))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(4.80000e+01))
        posep = findfirst(x->x=="C3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.10000e+02))
        ename = "E14"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eT1")
            arrset(ielftype,ie,iet_["eT1"])
        end
        vname = "X14"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(4.00000e+02))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(5.40000e+01))
        posep = findfirst(x->x=="C3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.00000e+02))
        ename = "E15"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eT1")
            arrset(ielftype,ie,iet_["eT1"])
        end
        vname = "X15"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(4.00000e+01))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(3.12200e+01))
        posep = findfirst(x->x=="C3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.30000e+02))
        ename = "E16"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eT1")
            arrset(ielftype,ie,iet_["eT1"])
        end
        vname = "X16"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(4.50000e+03))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(6.60000e+01))
        posep = findfirst(x->x=="C3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.22000e+02))
        ename = "E17"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eT1")
            arrset(ielftype,ie,iet_["eT1"])
        end
        vname = "X17"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(5.10000e+03))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(2.40000e+01))
        posep = findfirst(x->x=="C3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.10000e+02))
        ename = "E18"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eT1")
            arrset(ielftype,ie,iet_["eT1"])
        end
        vname = "X18"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(3.00000e+01))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(4.80000e+01))
        posep = findfirst(x->x=="C3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.18000e+02))
        ename = "E19"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eT4")
        arrset(ielftype,ie,iet_["eT4"])
        vname = "X19"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(4.84530e+02))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.12970e+02))
        ename = "E20"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eT1")
            arrset(ielftype,ie,iet_["eT1"])
        end
        vname = "X20"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(2.60000e+04))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(4.80000e+01))
        posep = findfirst(x->x=="C3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.10000e+02))
        ename = "E21"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eT4")
        arrset(ielftype,ie,iet_["eT4"])
        vname = "X21"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.86880e+02))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.60610e+02))
        ename = "E22"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eT1")
            arrset(ielftype,ie,iet_["eT1"])
        end
        vname = "X22"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(5.20000e+02))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(3.36100e+01))
        posep = findfirst(x->x=="C3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.30000e+02))
        ename = "E23"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eT1")
            arrset(ielftype,ie,iet_["eT1"])
        end
        vname = "X23"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(4.60000e+03))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(5.40000e+01))
        posep = findfirst(x->x=="C3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(9.50000e+01))
        ename = "E24"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eT1")
            arrset(ielftype,ie,iet_["eT1"])
        end
        vname = "X24"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(5.40000e+03))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(5.40000e+01))
        posep = findfirst(x->x=="C3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(9.50000e+01))
        ename = "E25"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eT1")
            arrset(ielftype,ie,iet_["eT1"])
        end
        vname = "X25"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(5.60000e+03))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.20000e+01))
        posep = findfirst(x->x=="C3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.10000e+02))
        ename = "E26"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eT1")
            arrset(ielftype,ie,iet_["eT1"])
        end
        vname = "X26"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(3.30000e+03))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.20000e+01))
        posep = findfirst(x->x=="C3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.10000e+02))
        ename = "E27"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eT1")
            arrset(ielftype,ie,iet_["eT1"])
        end
        vname = "X27"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(2.20000e+03))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(2.40000e+01))
        posep = findfirst(x->x=="C3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.24000e+02))
        ename = "E28"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eT1")
            arrset(ielftype,ie,iet_["eT1"])
        end
        vname = "X28"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.00000e+03))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(2.40000e+01))
        posep = findfirst(x->x=="C3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.10000e+02))
        ename = "E29"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eT1")
            arrset(ielftype,ie,iet_["eT1"])
        end
        vname = "X29"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(5.90000e+03))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(2.40000e+01))
        posep = findfirst(x->x=="C3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.13000e+02))
        ename = "E30"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eT1")
            arrset(ielftype,ie,iet_["eT1"])
        end
        vname = "X30"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(2.80000e+03))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(2.40000e+01))
        posep = findfirst(x->x=="C3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.13000e+02))
        ename = "E31"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eT1")
            arrset(ielftype,ie,iet_["eT1"])
        end
        vname = "X31"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(2.70000e+03))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.20000e+01))
        posep = findfirst(x->x=="C3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.10000e+02))
        ename = "E32"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eT1")
            arrset(ielftype,ie,iet_["eT1"])
        end
        vname = "X32"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(3.10000e+03))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(5.40000e+01))
        posep = findfirst(x->x=="C3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(9.50000e+01))
        ename = "E33"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eT1")
            arrset(ielftype,ie,iet_["eT1"])
        end
        vname = "X33"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.00000e+02))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.00000e+01))
        posep = findfirst(x->x=="C3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.00000e+02))
        ename = "E34"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eT1")
            arrset(ielftype,ie,iet_["eT1"])
        end
        vname = "X34"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(4.30000e+03))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(2.40000e+01))
        posep = findfirst(x->x=="C3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.13000e+02))
        ename = "E35"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eT1")
            arrset(ielftype,ie,iet_["eT1"])
        end
        vname = "X35"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(2.20000e+03))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(5.40000e+01))
        posep = findfirst(x->x=="C3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(9.50000e+01))
        ename = "E36"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eT1")
            arrset(ielftype,ie,iet_["eT1"])
        end
        vname = "X36"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.80000e+03))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(5.40000e+01))
        posep = findfirst(x->x=="C3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.10000e+02))
        ename = "E37"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eT1")
            arrset(ielftype,ie,iet_["eT1"])
        end
        vname = "X37"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.00000e+02))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(2.40000e+01))
        posep = findfirst(x->x=="C3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.10000e+02))
        ename = "E38"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eT1")
            arrset(ielftype,ie,iet_["eT1"])
        end
        vname = "X38"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.31000e+03))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(3.00000e+01))
        posep = findfirst(x->x=="C3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.00000e+02))
        ename = "E39"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eT1")
            arrset(ielftype,ie,iet_["eT1"])
        end
        vname = "X39"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(6.65000e+02))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(3.60000e+01))
        posep = findfirst(x->x=="C3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.20000e+02))
        ename = "E40"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eT1")
            arrset(ielftype,ie,iet_["eT1"])
        end
        vname = "X40"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.10000e+03))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(3.60000e+01))
        posep = findfirst(x->x=="C3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.20000e+02))
        ename = "E41"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eT1")
            arrset(ielftype,ie,iet_["eT1"])
        end
        vname = "X41"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(-2.00000e+02),Float64(2.00000e+02),Float64(-2.00000e+02)))
        posev = findfirst(x->x=="ARC",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="C1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(3.23000e+01))
        posep = findfirst(x->x=="C2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.00000e+01))
        posep = findfirst(x->x=="C3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.00000e+02))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["OBJ"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E2"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E3"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E4"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E5"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E6"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E7"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E8"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E9"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E10"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E11"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E12"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E13"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E14"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E15"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E16"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E17"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E18"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E19"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E20"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E21"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E22"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E23"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E24"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E25"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E26"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E27"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E28"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E29"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E30"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E31"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E32"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E33"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E34"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E35"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E36"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E37"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E38"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E39"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E40"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E41"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               -3.2393D+04
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
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-ONR2-MN-46-31"
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

    elseif action == "eT1"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        TMP = 850559.0e0/2.85*pbm.elpar[iel_][1]
        TMP = TMP/(pbm.elpar[iel_][3]^1.85)
        TMP = TMP/(pbm.elpar[iel_][2]^4.87)
        X = abs(EV_[1])
        XEXP = X^0.85
        f_   = TMP*X^2*XEXP
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 2.85*TMP*EV_[1]*XEXP
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 5.2725*TMP*XEXP
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eT4"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        EPS2 = 1.0e-14
        SQC1 = sqrt(pbm.elpar[iel_][1])
        X = min(EV_[1],SQC1)
        TMP = pbm.elpar[iel_][2]*(pbm.elpar[iel_][1]-X*X)
        TMP = sqrt(max(TMP,EPS2))
        TMP2 = sqrt(pbm.elpar[iel_][2])*asin(X/SQC1)
        f_   = 0.5*(-X*TMP-pbm.elpar[iel_][1]*TMP2)
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = -TMP
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = pbm.elpar[iel_][2]*X/TMP
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

