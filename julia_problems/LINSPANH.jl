function LINSPANH(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : LINSPANH
#    *********
# 
#    A linear network problem based on the spanish hydro-electric
#    reservoir management problem SPANHYD
# 
#    Source:
#    A partial specification of problem SPANHYD.
# 
#    SIF input: Ph. Toint, Sept 1990.
# 
#    classification = "C-LNR2-MN-97-33"
# 
#    Number of arcs = 97
#    Number of nodes
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "LINSPANH"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["NODES"] = 33
        v_["1"] = 1
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        for I = Int64(v_["1"]):Int64(v_["NODES"])
            ig,ig_,_ = s2mpj_ii("N"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"N"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        ngrp   = length(ig_)
        iv,ix_,_ = s2mpj_ii("X1",ix_)
        arrset(pb.xnames,iv,"X1")
        ig = ig_["OBJ"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("X1",ix_)
        arrset(pb.xnames,iv,"X1")
        ig = ig_["N32"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N1"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X2",ix_)
        arrset(pb.xnames,iv,"X2")
        ig = ig_["N32"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X3",ix_)
        arrset(pb.xnames,iv,"X3")
        ig = ig_["N32"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N3"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X4",ix_)
        arrset(pb.xnames,iv,"X4")
        ig = ig_["N32"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N4"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X5",ix_)
        arrset(pb.xnames,iv,"X5")
        ig = ig_["N32"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N5"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X6",ix_)
        arrset(pb.xnames,iv,"X6")
        ig = ig_["N32"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N6"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X7",ix_)
        arrset(pb.xnames,iv,"X7")
        ig = ig_["N32"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N7"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X8",ix_)
        arrset(pb.xnames,iv,"X8")
        ig = ig_["N32"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N8"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X9",ix_)
        arrset(pb.xnames,iv,"X9")
        ig = ig_["N32"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N9"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X10",ix_)
        arrset(pb.xnames,iv,"X10")
        ig = ig_["N32"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N10"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X11",ix_)
        arrset(pb.xnames,iv,"X11")
        ig = ig_["N1"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X12",ix_)
        arrset(pb.xnames,iv,"X12")
        ig = ig_["N2"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N3"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X13",ix_)
        arrset(pb.xnames,iv,"X13")
        ig = ig_["N3"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N5"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X14",ix_)
        arrset(pb.xnames,iv,"X14")
        ig = ig_["N4"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N5"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X15",ix_)
        arrset(pb.xnames,iv,"X15")
        ig = ig_["N5"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N9"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X16",ix_)
        arrset(pb.xnames,iv,"X16")
        ig = ig_["N6"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N7"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X17",ix_)
        arrset(pb.xnames,iv,"X17")
        ig = ig_["N7"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N8"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X18",ix_)
        arrset(pb.xnames,iv,"X18")
        ig = ig_["N8"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N9"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X19",ix_)
        arrset(pb.xnames,iv,"X19")
        ig = ig_["N9"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N10"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X20",ix_)
        arrset(pb.xnames,iv,"X20")
        ig = ig_["N10"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N31"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X21",ix_)
        arrset(pb.xnames,iv,"X21")
        ig = ig_["N1"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X22",ix_)
        arrset(pb.xnames,iv,"X22")
        ig = ig_["N2"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N3"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X23",ix_)
        arrset(pb.xnames,iv,"X23")
        ig = ig_["N3"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N9"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X24",ix_)
        arrset(pb.xnames,iv,"X24")
        ig = ig_["N4"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N9"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X25",ix_)
        arrset(pb.xnames,iv,"X25")
        ig = ig_["N6"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N7"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X26",ix_)
        arrset(pb.xnames,iv,"X26")
        ig = ig_["N7"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N8"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X27",ix_)
        arrset(pb.xnames,iv,"X27")
        ig = ig_["N8"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N9"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X28",ix_)
        arrset(pb.xnames,iv,"X28")
        ig = ig_["N9"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N10"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X29",ix_)
        arrset(pb.xnames,iv,"X29")
        ig = ig_["N10"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N31"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X30",ix_)
        arrset(pb.xnames,iv,"X30")
        ig = ig_["N1"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N11"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X31",ix_)
        arrset(pb.xnames,iv,"X31")
        ig = ig_["N2"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N12"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X32",ix_)
        arrset(pb.xnames,iv,"X32")
        ig = ig_["N3"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N13"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X33",ix_)
        arrset(pb.xnames,iv,"X33")
        ig = ig_["N4"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N14"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X34",ix_)
        arrset(pb.xnames,iv,"X34")
        ig = ig_["N5"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N15"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X35",ix_)
        arrset(pb.xnames,iv,"X35")
        ig = ig_["N6"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N16"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X36",ix_)
        arrset(pb.xnames,iv,"X36")
        ig = ig_["N7"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N17"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X37",ix_)
        arrset(pb.xnames,iv,"X37")
        ig = ig_["N8"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N18"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X38",ix_)
        arrset(pb.xnames,iv,"X38")
        ig = ig_["N9"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N19"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X39",ix_)
        arrset(pb.xnames,iv,"X39")
        ig = ig_["N10"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N20"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X40",ix_)
        arrset(pb.xnames,iv,"X40")
        ig = ig_["N11"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N12"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X41",ix_)
        arrset(pb.xnames,iv,"X41")
        ig = ig_["N12"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N13"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X42",ix_)
        arrset(pb.xnames,iv,"X42")
        ig = ig_["N13"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N15"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X43",ix_)
        arrset(pb.xnames,iv,"X43")
        ig = ig_["N14"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N15"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X44",ix_)
        arrset(pb.xnames,iv,"X44")
        ig = ig_["N15"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N19"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X45",ix_)
        arrset(pb.xnames,iv,"X45")
        ig = ig_["N16"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N17"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X46",ix_)
        arrset(pb.xnames,iv,"X46")
        ig = ig_["N17"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N18"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X47",ix_)
        arrset(pb.xnames,iv,"X47")
        ig = ig_["N18"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N19"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X48",ix_)
        arrset(pb.xnames,iv,"X48")
        ig = ig_["N19"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N20"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X49",ix_)
        arrset(pb.xnames,iv,"X49")
        ig = ig_["N20"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N31"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X50",ix_)
        arrset(pb.xnames,iv,"X50")
        ig = ig_["N11"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N12"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X51",ix_)
        arrset(pb.xnames,iv,"X51")
        ig = ig_["N12"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N13"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X52",ix_)
        arrset(pb.xnames,iv,"X52")
        ig = ig_["N13"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N19"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X53",ix_)
        arrset(pb.xnames,iv,"X53")
        ig = ig_["N14"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N19"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X54",ix_)
        arrset(pb.xnames,iv,"X54")
        ig = ig_["N16"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N17"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X55",ix_)
        arrset(pb.xnames,iv,"X55")
        ig = ig_["N17"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N18"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X56",ix_)
        arrset(pb.xnames,iv,"X56")
        ig = ig_["N18"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N19"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X57",ix_)
        arrset(pb.xnames,iv,"X57")
        ig = ig_["N19"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N20"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X58",ix_)
        arrset(pb.xnames,iv,"X58")
        ig = ig_["N20"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N31"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X59",ix_)
        arrset(pb.xnames,iv,"X59")
        ig = ig_["N11"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N21"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X60",ix_)
        arrset(pb.xnames,iv,"X60")
        ig = ig_["N12"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N22"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X61",ix_)
        arrset(pb.xnames,iv,"X61")
        ig = ig_["N13"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N23"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X62",ix_)
        arrset(pb.xnames,iv,"X62")
        ig = ig_["N14"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N24"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X63",ix_)
        arrset(pb.xnames,iv,"X63")
        ig = ig_["N15"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N25"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X64",ix_)
        arrset(pb.xnames,iv,"X64")
        ig = ig_["N16"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N26"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X65",ix_)
        arrset(pb.xnames,iv,"X65")
        ig = ig_["N17"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N27"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X66",ix_)
        arrset(pb.xnames,iv,"X66")
        ig = ig_["N18"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N28"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X67",ix_)
        arrset(pb.xnames,iv,"X67")
        ig = ig_["N19"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N29"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X68",ix_)
        arrset(pb.xnames,iv,"X68")
        ig = ig_["N20"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N30"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X69",ix_)
        arrset(pb.xnames,iv,"X69")
        ig = ig_["N21"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N22"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X70",ix_)
        arrset(pb.xnames,iv,"X70")
        ig = ig_["N22"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N23"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X71",ix_)
        arrset(pb.xnames,iv,"X71")
        ig = ig_["N23"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N25"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X72",ix_)
        arrset(pb.xnames,iv,"X72")
        ig = ig_["N24"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N25"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X73",ix_)
        arrset(pb.xnames,iv,"X73")
        ig = ig_["N25"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N29"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X74",ix_)
        arrset(pb.xnames,iv,"X74")
        ig = ig_["N26"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N27"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X75",ix_)
        arrset(pb.xnames,iv,"X75")
        ig = ig_["N27"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N28"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X76",ix_)
        arrset(pb.xnames,iv,"X76")
        ig = ig_["N28"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N29"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X77",ix_)
        arrset(pb.xnames,iv,"X77")
        ig = ig_["N29"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N30"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X78",ix_)
        arrset(pb.xnames,iv,"X78")
        ig = ig_["N30"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N31"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X79",ix_)
        arrset(pb.xnames,iv,"X79")
        ig = ig_["N21"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N22"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X80",ix_)
        arrset(pb.xnames,iv,"X80")
        ig = ig_["N22"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N23"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X81",ix_)
        arrset(pb.xnames,iv,"X81")
        ig = ig_["N23"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N29"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X82",ix_)
        arrset(pb.xnames,iv,"X82")
        ig = ig_["N24"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N29"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X83",ix_)
        arrset(pb.xnames,iv,"X83")
        ig = ig_["N26"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N27"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X84",ix_)
        arrset(pb.xnames,iv,"X84")
        ig = ig_["N27"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N28"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X85",ix_)
        arrset(pb.xnames,iv,"X85")
        ig = ig_["N28"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N29"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X86",ix_)
        arrset(pb.xnames,iv,"X86")
        ig = ig_["N29"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N30"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X87",ix_)
        arrset(pb.xnames,iv,"X87")
        ig = ig_["N30"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N31"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X88",ix_)
        arrset(pb.xnames,iv,"X88")
        ig = ig_["N21"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N33"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X89",ix_)
        arrset(pb.xnames,iv,"X89")
        ig = ig_["N22"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N33"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X90",ix_)
        arrset(pb.xnames,iv,"X90")
        ig = ig_["N23"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N33"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X91",ix_)
        arrset(pb.xnames,iv,"X91")
        ig = ig_["N24"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N33"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X92",ix_)
        arrset(pb.xnames,iv,"X92")
        ig = ig_["N25"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N33"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X93",ix_)
        arrset(pb.xnames,iv,"X93")
        ig = ig_["N26"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N33"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X94",ix_)
        arrset(pb.xnames,iv,"X94")
        ig = ig_["N27"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N33"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X95",ix_)
        arrset(pb.xnames,iv,"X95")
        ig = ig_["N28"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N33"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X96",ix_)
        arrset(pb.xnames,iv,"X96")
        ig = ig_["N29"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N33"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X97",ix_)
        arrset(pb.xnames,iv,"X97")
        ig = ig_["N30"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["N33"]
        pbm.A[ig,iv] += Float64(1.0)
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
        pbm.gconst[ig_["N1"]] = Float64(-5.13800e+01)
        pbm.gconst[ig_["N2"]] = Float64(-1.38400e+01)
        pbm.gconst[ig_["N3"]] = Float64(-2.58000)
        pbm.gconst[ig_["N4"]] = Float64(-2.19100e+01)
        pbm.gconst[ig_["N6"]] = Float64(-1.29700e+01)
        pbm.gconst[ig_["N8"]] = Float64(-2.89000)
        pbm.gconst[ig_["N9"]] = Float64(-2.08400e+01)
        pbm.gconst[ig_["N10"]] = Float64(-1.71400e+01)
        pbm.gconst[ig_["N11"]] = Float64(-3.20600e+01)
        pbm.gconst[ig_["N12"]] = Float64(-2.80000e-01)
        pbm.gconst[ig_["N13"]] = Float64(-4.20000)
        pbm.gconst[ig_["N14"]] = Float64(-4.83700e+01)
        pbm.gconst[ig_["N16"]] = Float64(-1.81300e+01)
        pbm.gconst[ig_["N18"]] = Float64(1.61000)
        pbm.gconst[ig_["N19"]] = Float64(-2.66000e+01)
        pbm.gconst[ig_["N20"]] = Float64(-1.87600e+01)
        pbm.gconst[ig_["N21"]] = Float64(-1.81300e+01)
        pbm.gconst[ig_["N24"]] = Float64(-1.81300e+01)
        pbm.gconst[ig_["N26"]] = Float64(-9.10000)
        pbm.gconst[ig_["N28"]] = Float64(5.81000)
        pbm.gconst[ig_["N29"]] = Float64(-9.10000)
        pbm.gconst[ig_["N30"]] = Float64(-6.02000)
        pbm.gconst[ig_["N31"]] = Float64(6.08350e+02)
        pbm.gconst[ig_["N32"]] = Float64(-4.62634e+03)
        pbm.gconst[ig_["N33"]] = Float64(4.36300e+03)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xupper = fill(3.02400e+03,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        pb.xlower[ix_["X1"]] = 7.70000e+01
        pb.xupper[ix_["X1"]] = 7.70100e+01
        pb.xlower[ix_["X2"]] = 1.12452e+03
        pb.xupper[ix_["X2"]] = 1.12453e+03
        pb.xlower[ix_["X3"]] = 1.58000e+02
        pb.xupper[ix_["X3"]] = 1.58010e+02
        pb.xlower[ix_["X4"]] = 1.60000e+01
        pb.xupper[ix_["X4"]] = 1.60100e+01
        pb.xlower[ix_["X5"]] = 0.00000
        pb.xupper[ix_["X5"]] = 0.00000
        pb.xlower[ix_["X6"]] = 7.83650e+02
        pb.xupper[ix_["X6"]] = 7.83660e+02
        pb.xlower[ix_["X7"]] = 1.10000e+01
        pb.xupper[ix_["X7"]] = 1.10100e+01
        pb.xlower[ix_["X8"]] = 4.90000e+01
        pb.xupper[ix_["X8"]] = 4.90100e+01
        pb.xlower[ix_["X9"]] = 2.15517e+03
        pb.xupper[ix_["X9"]] = 2.15518e+03
        pb.xlower[ix_["X10"]] = 2.52000e+02
        pb.xupper[ix_["X10"]] = 2.52010e+02
        pb.xupper[ix_["X11"]] = 3.97840e+02
        pb.xupper[ix_["X12"]] = 2.22320e+02
        pb.xupper[ix_["X13"]] = 2.05630e+02
        pb.xupper[ix_["X14"]] = 2.05630e+02
        pb.xupper[ix_["X15"]] = 2.05630e+02
        pb.xupper[ix_["X16"]] = 1.24830e+02
        pb.xupper[ix_["X17"]] = 1.27010e+02
        pb.xupper[ix_["X18"]] = 6.10800e+01
        pb.xupper[ix_["X19"]] = 6.14840e+02
        pb.xupper[ix_["X20"]] = 7.78080e+02
        pb.xupper[ix_["X25"]] = 7.25760e+03
        pb.xupper[ix_["X26"]] = 1.20960e+03
        pb.xupper[ix_["X27"]] = 9.07200e+02
        pb.xupper[ix_["X28"]] = 7.25760e+03
        pb.xupper[ix_["X29"]] = 7.25760e+03
        pb.xlower[ix_["X30"]] = 7.70000e+01
        pb.xupper[ix_["X30"]] = 7.70000e+01
        pb.xlower[ix_["X31"]] = 4.03400e+02
        pb.xupper[ix_["X31"]] = 1.31200e+03
        pb.xlower[ix_["X32"]] = 1.58000e+02
        pb.xupper[ix_["X32"]] = 1.58000e+02
        pb.xlower[ix_["X33"]] = 1.60000e+01
        pb.xupper[ix_["X33"]] = 1.60000e+01
        pb.xlower[ix_["X34"]] = 0.00000
        pb.xupper[ix_["X34"]] = 0.00000
        pb.xlower[ix_["X35"]] = 5.02000e+02
        pb.xupper[ix_["X35"]] = 9.28460e+02
        pb.xlower[ix_["X36"]] = 1.10000e+01
        pb.xupper[ix_["X36"]] = 1.10000e+01
        pb.xlower[ix_["X37"]] = 4.90000e+01
        pb.xupper[ix_["X37"]] = 4.90000e+01
        pb.xlower[ix_["X38"]] = 9.15300e+02
        pb.xupper[ix_["X38"]] = 2.61160e+03
        pb.xlower[ix_["X39"]] = 2.52000e+02
        pb.xupper[ix_["X39"]] = 2.52000e+02
        pb.xupper[ix_["X40"]] = 3.97840e+02
        pb.xupper[ix_["X41"]] = 2.22320e+02
        pb.xupper[ix_["X42"]] = 2.05630e+02
        pb.xupper[ix_["X43"]] = 2.05630e+02
        pb.xupper[ix_["X44"]] = 2.05630e+02
        pb.xupper[ix_["X45"]] = 1.24830e+02
        pb.xupper[ix_["X46"]] = 1.27010e+02
        pb.xupper[ix_["X47"]] = 6.10800e+01
        pb.xupper[ix_["X48"]] = 6.14840e+02
        pb.xupper[ix_["X49"]] = 7.78080e+02
        pb.xupper[ix_["X54"]] = 7.25760e+03
        pb.xupper[ix_["X55"]] = 1.20960e+03
        pb.xupper[ix_["X56"]] = 9.07200e+02
        pb.xupper[ix_["X57"]] = 7.25760e+03
        pb.xupper[ix_["X58"]] = 7.25760e+03
        pb.xlower[ix_["X59"]] = 7.70000e+01
        pb.xupper[ix_["X59"]] = 7.70000e+01
        pb.xlower[ix_["X60"]] = 4.03400e+02
        pb.xupper[ix_["X60"]] = 1.31200e+03
        pb.xlower[ix_["X61"]] = 1.58000e+02
        pb.xupper[ix_["X61"]] = 1.58000e+02
        pb.xlower[ix_["X62"]] = 1.60000e+01
        pb.xupper[ix_["X62"]] = 1.60000e+01
        pb.xlower[ix_["X63"]] = 0.00000
        pb.xupper[ix_["X63"]] = 0.00000
        pb.xlower[ix_["X64"]] = 5.05640e+02
        pb.xupper[ix_["X64"]] = 9.28460e+02
        pb.xlower[ix_["X65"]] = 1.10000e+01
        pb.xupper[ix_["X65"]] = 1.10000e+01
        pb.xlower[ix_["X66"]] = 4.90000e+01
        pb.xupper[ix_["X66"]] = 4.90000e+01
        pb.xlower[ix_["X67"]] = 9.15300e+02
        pb.xupper[ix_["X67"]] = 2.61160e+03
        pb.xlower[ix_["X68"]] = 2.52000e+02
        pb.xupper[ix_["X68"]] = 2.52000e+02
        pb.xupper[ix_["X69"]] = 3.97840e+02
        pb.xupper[ix_["X70"]] = 2.22320e+02
        pb.xupper[ix_["X71"]] = 2.05630e+02
        pb.xupper[ix_["X72"]] = 2.05630e+02
        pb.xupper[ix_["X73"]] = 2.05630e+02
        pb.xupper[ix_["X74"]] = 1.24830e+02
        pb.xupper[ix_["X75"]] = 1.27010e+02
        pb.xupper[ix_["X76"]] = 6.10800e+01
        pb.xupper[ix_["X77"]] = 6.14840e+02
        pb.xupper[ix_["X78"]] = 7.78080e+02
        pb.xupper[ix_["X83"]] = 7.25760e+03
        pb.xupper[ix_["X84"]] = 1.20960e+03
        pb.xupper[ix_["X85"]] = 9.07200e+02
        pb.xupper[ix_["X86"]] = 7.25760e+03
        pb.xupper[ix_["X87"]] = 7.25760e+03
        pb.xlower[ix_["X88"]] = 7.70000e+01
        pb.xupper[ix_["X88"]] = 7.70100e+01
        pb.xlower[ix_["X89"]] = 1.10000e+03
        pb.xupper[ix_["X89"]] = 1.10001e+03
        pb.xlower[ix_["X90"]] = 1.58000e+02
        pb.xupper[ix_["X90"]] = 1.58010e+02
        pb.xlower[ix_["X91"]] = 1.60000e+01
        pb.xupper[ix_["X91"]] = 1.60100e+01
        pb.xlower[ix_["X92"]] = 0.00000
        pb.xupper[ix_["X92"]] = 0.00000
        pb.xlower[ix_["X93"]] = 7.00000e+02
        pb.xupper[ix_["X93"]] = 7.00010e+02
        pb.xlower[ix_["X94"]] = 1.10000e+01
        pb.xupper[ix_["X94"]] = 1.10100e+01
        pb.xlower[ix_["X95"]] = 4.90000e+01
        pb.xupper[ix_["X95"]] = 4.90100e+01
        pb.xlower[ix_["X96"]] = 2.00000e+03
        pb.xupper[ix_["X96"]] = 2.00001e+03
        pb.xlower[ix_["X97"]] = 2.52000e+02
        pb.xupper[ix_["X97"]] = 2.52010e+02
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(0.0),pb.n)
        if haskey(ix_,"X1")
            pb.x0[ix_["X1"]] = Float64(7.70000e+01)
        else
            pb.y0[findfirst(x->x==ig_["X1"],pbm.congrps)] = Float64(7.70000e+01)
        end
        if haskey(ix_,"X2")
            pb.x0[ix_["X2"]] = Float64(1.12452e+03)
        else
            pb.y0[findfirst(x->x==ig_["X2"],pbm.congrps)] = Float64(1.12452e+03)
        end
        if haskey(ix_,"X3")
            pb.x0[ix_["X3"]] = Float64(1.58000e+02)
        else
            pb.y0[findfirst(x->x==ig_["X3"],pbm.congrps)] = Float64(1.58000e+02)
        end
        if haskey(ix_,"X4")
            pb.x0[ix_["X4"]] = Float64(1.60000e+01)
        else
            pb.y0[findfirst(x->x==ig_["X4"],pbm.congrps)] = Float64(1.60000e+01)
        end
        if haskey(ix_,"X6")
            pb.x0[ix_["X6"]] = Float64(7.83650e+02)
        else
            pb.y0[findfirst(x->x==ig_["X6"],pbm.congrps)] = Float64(7.83650e+02)
        end
        if haskey(ix_,"X7")
            pb.x0[ix_["X7"]] = Float64(1.10000e+01)
        else
            pb.y0[findfirst(x->x==ig_["X7"],pbm.congrps)] = Float64(1.10000e+01)
        end
        if haskey(ix_,"X8")
            pb.x0[ix_["X8"]] = Float64(4.90000e+01)
        else
            pb.y0[findfirst(x->x==ig_["X8"],pbm.congrps)] = Float64(4.90000e+01)
        end
        if haskey(ix_,"X9")
            pb.x0[ix_["X9"]] = Float64(2.15517e+03)
        else
            pb.y0[findfirst(x->x==ig_["X9"],pbm.congrps)] = Float64(2.15517e+03)
        end
        if haskey(ix_,"X10")
            pb.x0[ix_["X10"]] = Float64(2.52000e+02)
        else
            pb.y0[findfirst(x->x==ig_["X10"],pbm.congrps)] = Float64(2.52000e+02)
        end
        if haskey(ix_,"X11")
            pb.x0[ix_["X11"]] = Float64(5.13800e+01)
        else
            pb.y0[findfirst(x->x==ig_["X11"],pbm.congrps)] = Float64(5.13800e+01)
        end
        if haskey(ix_,"X12")
            pb.x0[ix_["X12"]] = Float64(1.40210e+02)
        else
            pb.y0[findfirst(x->x==ig_["X12"],pbm.congrps)] = Float64(1.40210e+02)
        end
        if haskey(ix_,"X13")
            pb.x0[ix_["X13"]] = Float64(1.42790e+02)
        else
            pb.y0[findfirst(x->x==ig_["X13"],pbm.congrps)] = Float64(1.42790e+02)
        end
        if haskey(ix_,"X14")
            pb.x0[ix_["X14"]] = Float64(2.19100e+01)
        else
            pb.y0[findfirst(x->x==ig_["X14"],pbm.congrps)] = Float64(2.19100e+01)
        end
        if haskey(ix_,"X15")
            pb.x0[ix_["X15"]] = Float64(1.64700e+02)
        else
            pb.y0[findfirst(x->x==ig_["X15"],pbm.congrps)] = Float64(1.64700e+02)
        end
        if haskey(ix_,"X16")
            pb.x0[ix_["X16"]] = Float64(5.81900e+01)
        else
            pb.y0[findfirst(x->x==ig_["X16"],pbm.congrps)] = Float64(5.81900e+01)
        end
        if haskey(ix_,"X17")
            pb.x0[ix_["X17"]] = Float64(5.81900e+01)
        else
            pb.y0[findfirst(x->x==ig_["X17"],pbm.congrps)] = Float64(5.81900e+01)
        end
        if haskey(ix_,"X18")
            pb.x0[ix_["X18"]] = Float64(6.10800e+01)
        else
            pb.y0[findfirst(x->x==ig_["X18"],pbm.congrps)] = Float64(6.10800e+01)
        end
        if haskey(ix_,"X19")
            pb.x0[ix_["X19"]] = Float64(5.66430e+02)
        else
            pb.y0[findfirst(x->x==ig_["X19"],pbm.congrps)] = Float64(5.66430e+02)
        end
        if haskey(ix_,"X20")
            pb.x0[ix_["X20"]] = Float64(5.83570e+02)
        else
            pb.y0[findfirst(x->x==ig_["X20"],pbm.congrps)] = Float64(5.83570e+02)
        end
        if haskey(ix_,"X30")
            pb.x0[ix_["X30"]] = Float64(7.70000e+01)
        else
            pb.y0[findfirst(x->x==ig_["X30"],pbm.congrps)] = Float64(7.70000e+01)
        end
        if haskey(ix_,"X31")
            pb.x0[ix_["X31"]] = Float64(1.04953e+03)
        else
            pb.y0[findfirst(x->x==ig_["X31"],pbm.congrps)] = Float64(1.04953e+03)
        end
        if haskey(ix_,"X32")
            pb.x0[ix_["X32"]] = Float64(1.58000e+02)
        else
            pb.y0[findfirst(x->x==ig_["X32"],pbm.congrps)] = Float64(1.58000e+02)
        end
        if haskey(ix_,"X33")
            pb.x0[ix_["X33"]] = Float64(1.60000e+01)
        else
            pb.y0[findfirst(x->x==ig_["X33"],pbm.congrps)] = Float64(1.60000e+01)
        end
        if haskey(ix_,"X35")
            pb.x0[ix_["X35"]] = Float64(7.38430e+02)
        else
            pb.y0[findfirst(x->x==ig_["X35"],pbm.congrps)] = Float64(7.38430e+02)
        end
        if haskey(ix_,"X36")
            pb.x0[ix_["X36"]] = Float64(1.10000e+01)
        else
            pb.y0[findfirst(x->x==ig_["X36"],pbm.congrps)] = Float64(1.10000e+01)
        end
        if haskey(ix_,"X37")
            pb.x0[ix_["X37"]] = Float64(4.90000e+01)
        else
            pb.y0[findfirst(x->x==ig_["X37"],pbm.congrps)] = Float64(4.90000e+01)
        end
        if haskey(ix_,"X38")
            pb.x0[ix_["X38"]] = Float64(1.83536e+03)
        else
            pb.y0[findfirst(x->x==ig_["X38"],pbm.congrps)] = Float64(1.83536e+03)
        end
        if haskey(ix_,"X39")
            pb.x0[ix_["X39"]] = Float64(2.52000e+02)
        else
            pb.y0[findfirst(x->x==ig_["X39"],pbm.congrps)] = Float64(2.52000e+02)
        end
        if haskey(ix_,"X40")
            pb.x0[ix_["X40"]] = Float64(3.20600e+01)
        else
            pb.y0[findfirst(x->x==ig_["X40"],pbm.congrps)] = Float64(3.20600e+01)
        end
        if haskey(ix_,"X42")
            pb.x0[ix_["X42"]] = Float64(4.20000)
        else
            pb.y0[findfirst(x->x==ig_["X42"],pbm.congrps)] = Float64(4.20000)
        end
        if haskey(ix_,"X43")
            pb.x0[ix_["X43"]] = Float64(4.83700e+01)
        else
            pb.y0[findfirst(x->x==ig_["X43"],pbm.congrps)] = Float64(4.83700e+01)
        end
        if haskey(ix_,"X44")
            pb.x0[ix_["X44"]] = Float64(5.25700e+01)
        else
            pb.y0[findfirst(x->x==ig_["X44"],pbm.congrps)] = Float64(5.25700e+01)
        end
        if haskey(ix_,"X45")
            pb.x0[ix_["X45"]] = Float64(5.98500e+01)
        else
            pb.y0[findfirst(x->x==ig_["X45"],pbm.congrps)] = Float64(5.98500e+01)
        end
        if haskey(ix_,"X46")
            pb.x0[ix_["X46"]] = Float64(5.98500e+01)
        else
            pb.y0[findfirst(x->x==ig_["X46"],pbm.congrps)] = Float64(5.98500e+01)
        end
        if haskey(ix_,"X47")
            pb.x0[ix_["X47"]] = Float64(5.82400e+01)
        else
            pb.y0[findfirst(x->x==ig_["X47"],pbm.congrps)] = Float64(5.82400e+01)
        end
        if haskey(ix_,"X49")
            pb.x0[ix_["X49"]] = Float64(1.87600e+01)
        else
            pb.y0[findfirst(x->x==ig_["X49"],pbm.congrps)] = Float64(1.87600e+01)
        end
        if haskey(ix_,"X59")
            pb.x0[ix_["X59"]] = Float64(7.70000e+01)
        else
            pb.y0[findfirst(x->x==ig_["X59"],pbm.congrps)] = Float64(7.70000e+01)
        end
        if haskey(ix_,"X60")
            pb.x0[ix_["X60"]] = Float64(1.08187e+03)
        else
            pb.y0[findfirst(x->x==ig_["X60"],pbm.congrps)] = Float64(1.08187e+03)
        end
        if haskey(ix_,"X61")
            pb.x0[ix_["X61"]] = Float64(1.58000e+02)
        else
            pb.y0[findfirst(x->x==ig_["X61"],pbm.congrps)] = Float64(1.58000e+02)
        end
        if haskey(ix_,"X62")
            pb.x0[ix_["X62"]] = Float64(1.60000e+01)
        else
            pb.y0[findfirst(x->x==ig_["X62"],pbm.congrps)] = Float64(1.60000e+01)
        end
        if haskey(ix_,"X64")
            pb.x0[ix_["X64"]] = Float64(6.96710e+02)
        else
            pb.y0[findfirst(x->x==ig_["X64"],pbm.congrps)] = Float64(6.96710e+02)
        end
        if haskey(ix_,"X65")
            pb.x0[ix_["X65"]] = Float64(1.10000e+01)
        else
            pb.y0[findfirst(x->x==ig_["X65"],pbm.congrps)] = Float64(1.10000e+01)
        end
        if haskey(ix_,"X66")
            pb.x0[ix_["X66"]] = Float64(4.90000e+01)
        else
            pb.y0[findfirst(x->x==ig_["X66"],pbm.congrps)] = Float64(4.90000e+01)
        end
        if haskey(ix_,"X67")
            pb.x0[ix_["X67"]] = Float64(1.97277e+03)
        else
            pb.y0[findfirst(x->x==ig_["X67"],pbm.congrps)] = Float64(1.97277e+03)
        end
        if haskey(ix_,"X68")
            pb.x0[ix_["X68"]] = Float64(2.52000e+02)
        else
            pb.y0[findfirst(x->x==ig_["X68"],pbm.congrps)] = Float64(2.52000e+02)
        end
        if haskey(ix_,"X69")
            pb.x0[ix_["X69"]] = Float64(1.81300e+01)
        else
            pb.y0[findfirst(x->x==ig_["X69"],pbm.congrps)] = Float64(1.81300e+01)
        end
        if haskey(ix_,"X72")
            pb.x0[ix_["X72"]] = Float64(1.81300e+01)
        else
            pb.y0[findfirst(x->x==ig_["X72"],pbm.congrps)] = Float64(1.81300e+01)
        end
        if haskey(ix_,"X73")
            pb.x0[ix_["X73"]] = Float64(1.81300e+01)
        else
            pb.y0[findfirst(x->x==ig_["X73"],pbm.congrps)] = Float64(1.81300e+01)
        end
        if haskey(ix_,"X74")
            pb.x0[ix_["X74"]] = Float64(5.81000)
        else
            pb.y0[findfirst(x->x==ig_["X74"],pbm.congrps)] = Float64(5.81000)
        end
        if haskey(ix_,"X75")
            pb.x0[ix_["X75"]] = Float64(5.81000)
        else
            pb.y0[findfirst(x->x==ig_["X75"],pbm.congrps)] = Float64(5.81000)
        end
        if haskey(ix_,"X78")
            pb.x0[ix_["X78"]] = Float64(6.02000)
        else
            pb.y0[findfirst(x->x==ig_["X78"],pbm.congrps)] = Float64(6.02000)
        end
        if haskey(ix_,"X88")
            pb.x0[ix_["X88"]] = Float64(7.70000e+01)
        else
            pb.y0[findfirst(x->x==ig_["X88"],pbm.congrps)] = Float64(7.70000e+01)
        end
        if haskey(ix_,"X89")
            pb.x0[ix_["X89"]] = Float64(1.10000e+03)
        else
            pb.y0[findfirst(x->x==ig_["X89"],pbm.congrps)] = Float64(1.10000e+03)
        end
        if haskey(ix_,"X90")
            pb.x0[ix_["X90"]] = Float64(1.58000e+02)
        else
            pb.y0[findfirst(x->x==ig_["X90"],pbm.congrps)] = Float64(1.58000e+02)
        end
        if haskey(ix_,"X91")
            pb.x0[ix_["X91"]] = Float64(1.60000e+01)
        else
            pb.y0[findfirst(x->x==ig_["X91"],pbm.congrps)] = Float64(1.60000e+01)
        end
        if haskey(ix_,"X93")
            pb.x0[ix_["X93"]] = Float64(7.00000e+02)
        else
            pb.y0[findfirst(x->x==ig_["X93"],pbm.congrps)] = Float64(7.00000e+02)
        end
        if haskey(ix_,"X94")
            pb.x0[ix_["X94"]] = Float64(1.10000e+01)
        else
            pb.y0[findfirst(x->x==ig_["X94"],pbm.congrps)] = Float64(1.10000e+01)
        end
        if haskey(ix_,"X95")
            pb.x0[ix_["X95"]] = Float64(4.90000e+01)
        else
            pb.y0[findfirst(x->x==ig_["X95"],pbm.congrps)] = Float64(4.90000e+01)
        end
        if haskey(ix_,"X96")
            pb.x0[ix_["X96"]] = Float64(2.00000e+03)
        else
            pb.y0[findfirst(x->x==ig_["X96"],pbm.congrps)] = Float64(2.00000e+03)
        end
        if haskey(ix_,"X97")
            pb.x0[ix_["X97"]] = Float64(2.52000e+02)
        else
            pb.y0[findfirst(x->x==ig_["X97"],pbm.congrps)] = Float64(2.52000e+02)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN                77.0
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
        pb.pbclass = "C-LNR2-MN-97-33"
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

