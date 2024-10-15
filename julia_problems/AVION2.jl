function AVION2(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : AVION2
#    *********
# 
#    Dassault France avion (airplane design) problem
# 
#    SIF input:  A. R. Conn, June 1993.
# 
#    classification = "C-OLR2-RN-49-15"
# 
#    Define useful parameters
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 6 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "AVION2"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["790"] = 790.0
        v_["1/790"] = 1.0/v_["790"]
        v_["24000"] = 24000.0
        v_["1/24000"] = 1.0/v_["24000"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        iv,ix_,_ = s2mpj_ii("SR",ix_)
        arrset(pb.xnames,iv,"SR")
        iv,ix_,_ = s2mpj_ii("LR",ix_)
        arrset(pb.xnames,iv,"LR")
        iv,ix_,_ = s2mpj_ii("PK",ix_)
        arrset(pb.xnames,iv,"PK")
        iv,ix_,_ = s2mpj_ii("EF",ix_)
        arrset(pb.xnames,iv,"EF")
        iv,ix_,_ = s2mpj_ii("SX",ix_)
        arrset(pb.xnames,iv,"SX")
        iv,ix_,_ = s2mpj_ii("LX",ix_)
        arrset(pb.xnames,iv,"LX")
        iv,ix_,_ = s2mpj_ii("SD",ix_)
        arrset(pb.xnames,iv,"SD")
        iv,ix_,_ = s2mpj_ii("SK",ix_)
        arrset(pb.xnames,iv,"SK")
        iv,ix_,_ = s2mpj_ii("ST",ix_)
        arrset(pb.xnames,iv,"ST")
        iv,ix_,_ = s2mpj_ii("SF",ix_)
        arrset(pb.xnames,iv,"SF")
        iv,ix_,_ = s2mpj_ii("LF",ix_)
        arrset(pb.xnames,iv,"LF")
        iv,ix_,_ = s2mpj_ii("AM",ix_)
        arrset(pb.xnames,iv,"AM")
        iv,ix_,_ = s2mpj_ii("CA",ix_)
        arrset(pb.xnames,iv,"CA")
        iv,ix_,_ = s2mpj_ii("CB",ix_)
        arrset(pb.xnames,iv,"CB")
        iv,ix_,_ = s2mpj_ii("SO",ix_)
        arrset(pb.xnames,iv,"SO")
        iv,ix_,_ = s2mpj_ii("SS",ix_)
        arrset(pb.xnames,iv,"SS")
        iv,ix_,_ = s2mpj_ii("IMPDER",ix_)
        arrset(pb.xnames,iv,"IMPDER")
        iv,ix_,_ = s2mpj_ii("IMPK",ix_)
        arrset(pb.xnames,iv,"IMPK")
        iv,ix_,_ = s2mpj_ii("IMPFUS",ix_)
        arrset(pb.xnames,iv,"IMPFUS")
        iv,ix_,_ = s2mpj_ii("QI",ix_)
        arrset(pb.xnames,iv,"QI")
        iv,ix_,_ = s2mpj_ii("PT",ix_)
        arrset(pb.xnames,iv,"PT")
        iv,ix_,_ = s2mpj_ii("MV",ix_)
        arrset(pb.xnames,iv,"MV")
        iv,ix_,_ = s2mpj_ii("MC",ix_)
        arrset(pb.xnames,iv,"MC")
        iv,ix_,_ = s2mpj_ii("MD",ix_)
        arrset(pb.xnames,iv,"MD")
        iv,ix_,_ = s2mpj_ii("PD",ix_)
        arrset(pb.xnames,iv,"PD")
        iv,ix_,_ = s2mpj_ii("NS",ix_)
        arrset(pb.xnames,iv,"NS")
        iv,ix_,_ = s2mpj_ii("VS",ix_)
        arrset(pb.xnames,iv,"VS")
        iv,ix_,_ = s2mpj_ii("CR",ix_)
        arrset(pb.xnames,iv,"CR")
        iv,ix_,_ = s2mpj_ii("PM",ix_)
        arrset(pb.xnames,iv,"PM")
        iv,ix_,_ = s2mpj_ii("DV",ix_)
        arrset(pb.xnames,iv,"DV")
        iv,ix_,_ = s2mpj_ii("MZ",ix_)
        arrset(pb.xnames,iv,"MZ")
        iv,ix_,_ = s2mpj_ii("VN",ix_)
        arrset(pb.xnames,iv,"VN")
        iv,ix_,_ = s2mpj_ii("QV",ix_)
        arrset(pb.xnames,iv,"QV")
        iv,ix_,_ = s2mpj_ii("QF",ix_)
        arrset(pb.xnames,iv,"QF")
        iv,ix_,_ = s2mpj_ii("IMPTRAIN",ix_)
        arrset(pb.xnames,iv,"IMPTRAIN")
        iv,ix_,_ = s2mpj_ii("IMPMOT",ix_)
        arrset(pb.xnames,iv,"IMPMOT")
        iv,ix_,_ = s2mpj_ii("IMPNMOT",ix_)
        arrset(pb.xnames,iv,"IMPNMOT")
        iv,ix_,_ = s2mpj_ii("IMPPET",ix_)
        arrset(pb.xnames,iv,"IMPPET")
        iv,ix_,_ = s2mpj_ii("IMPPIL",ix_)
        arrset(pb.xnames,iv,"IMPPIL")
        iv,ix_,_ = s2mpj_ii("IMPCAN",ix_)
        arrset(pb.xnames,iv,"IMPCAN")
        iv,ix_,_ = s2mpj_ii("IMPSNA",ix_)
        arrset(pb.xnames,iv,"IMPSNA")
        iv,ix_,_ = s2mpj_ii("MS",ix_)
        arrset(pb.xnames,iv,"MS")
        iv,ix_,_ = s2mpj_ii("EL",ix_)
        arrset(pb.xnames,iv,"EL")
        iv,ix_,_ = s2mpj_ii("DE",ix_)
        arrset(pb.xnames,iv,"DE")
        iv,ix_,_ = s2mpj_ii("DS",ix_)
        arrset(pb.xnames,iv,"DS")
        iv,ix_,_ = s2mpj_ii("IMPVOIL",ix_)
        arrset(pb.xnames,iv,"IMPVOIL")
        iv,ix_,_ = s2mpj_ii("NM",ix_)
        arrset(pb.xnames,iv,"NM")
        iv,ix_,_ = s2mpj_ii("NP",ix_)
        arrset(pb.xnames,iv,"NP")
        iv,ix_,_ = s2mpj_ii("NG",ix_)
        arrset(pb.xnames,iv,"NG")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("E1",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"E1")
        iv = ix_["SD"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["SR"]
        pbm.A[ig,iv] += Float64(-0.13)
        ig,ig_,_ = s2mpj_ii("E2",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"E2")
        iv = ix_["SX"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["SR"]
        pbm.A[ig,iv] += Float64(-0.7)
        ig,ig_,_ = s2mpj_ii("E3",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"E3")
        iv = ix_["LX"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["LR"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("E4",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["SK"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("E5",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"E5")
        iv = ix_["SF"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["ST"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["SD"]
        pbm.A[ig,iv] += Float64(-2.0)
        iv = ix_["SX"]
        pbm.A[ig,iv] += Float64(-2.0)
        iv = ix_["SK"]
        pbm.A[ig,iv] += Float64(-2.0)
        ig,ig_,_ = s2mpj_ii("E6",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["CA"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("E7",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["AM"]
        pbm.A[ig,iv] += Float64(-2.0)
        iv = ix_["SO"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["SS"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("E8",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["AM"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("E9",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["IMPDER"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["SD"]
        pbm.A[ig,iv] += Float64(-27.5)
        ig,ig_,_ = s2mpj_ii("E10",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["IMPK"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["SK"]
        pbm.A[ig,iv] += Float64(-70.0)
        ig,ig_,_ = s2mpj_ii("E11",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"E11")
        iv = ix_["IMPFUS"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["SF"]
        pbm.A[ig,iv] += Float64(-20.0)
        ig,ig_,_ = s2mpj_ii("E12",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"E12")
        iv = ix_["MD"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["MV"]
        pbm.A[ig,iv] += Float64(-2.0)
        ig,ig_,_ = s2mpj_ii("E13",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["QI"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("E14",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["PT"]
        pbm.A[ig,iv] += Float64(1000.0)
        ig,ig_,_ = s2mpj_ii("E15",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"E15")
        iv = ix_["QF"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["QI"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["QV"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("E16",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["VN"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["VS"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["QF"]
        pbm.A[ig,iv] += Float64(v_["1/790"])
        ig,ig_,_ = s2mpj_ii("E17",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"E17")
        iv = ix_["IMPTRAIN"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["MV"]
        pbm.A[ig,iv] += Float64(-0.137)
        ig,ig_,_ = s2mpj_ii("E18",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["IMPMOT"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("E19",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"E19")
        iv = ix_["IMPNMOT"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["NM"]
        pbm.A[ig,iv] += Float64(-35.0)
        ig,ig_,_ = s2mpj_ii("E20",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"E20")
        iv = ix_["IMPPET"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["QI"]
        pbm.A[ig,iv] += Float64(-0.043)
        ig,ig_,_ = s2mpj_ii("E21",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"E21")
        iv = ix_["IMPPIL"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["NP"]
        pbm.A[ig,iv] += Float64(-200.0)
        ig,ig_,_ = s2mpj_ii("E22",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"E22")
        iv = ix_["IMPCAN"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["NG"]
        pbm.A[ig,iv] += Float64(-120.0)
        ig,ig_,_ = s2mpj_ii("E23",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"E23")
        iv = ix_["IMPSNA"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["NS"]
        pbm.A[ig,iv] += Float64(-300.0)
        ig,ig_,_ = s2mpj_ii("E24",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"E24")
        iv = ix_["MC"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["MV"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["NP"]
        pbm.A[ig,iv] += Float64(95.0)
        iv = ix_["NG"]
        pbm.A[ig,iv] += Float64(70.0)
        iv = ix_["NM"]
        pbm.A[ig,iv] += Float64(660.0)
        iv = ix_["QI"]
        pbm.A[ig,iv] += Float64(0.5)
        ig,ig_,_ = s2mpj_ii("E25",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"E25")
        iv = ix_["MZ"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["IMPTRAIN"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["IMPNMOT"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["IMPPET"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["IMPPIL"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["IMPCAN"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["IMPSNA"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("E26",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["ST"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("E27",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["SR"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("E28",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["QV"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("E29",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["SO"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("E30",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["SS"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("E31",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["CB"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("E32",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["IMPVOIL"]
        pbm.A[ig,iv] += Float64(1.0)
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
        pbm.gconst[ig_["E13"]] = Float64(1000.0)
        pbm.gconst[ig_["E16"]] = Float64(-2.0)
        pbm.gconst[ig_["E23"]] = Float64(400.0)
        pbm.gconst[ig_["E24"]] = Float64(380.0)
        pbm.gconst[ig_["E25"]] = Float64(-290.0)
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower[ix_["SR"]] = 10.0
        pb.xlower[ix_["LR"]] = 0.0
        pb.xlower[ix_["PK"]] = 0.0
        pb.xlower[ix_["EF"]] = 0.0
        pb.xlower[ix_["SX"]] = 7.0
        pb.xlower[ix_["LX"]] = 1.5
        pb.xlower[ix_["SD"]] = 2.0
        pb.xlower[ix_["SK"]] = 2.0
        pb.xlower[ix_["ST"]] = 30.0
        pb.xlower[ix_["SF"]] = 20.0
        pb.xlower[ix_["LF"]] = 0.001
        pb.xlower[ix_["AM"]] = 0.0
        pb.xlower[ix_["CA"]] = -0.2
        pb.xlower[ix_["CB"]] = 0.1
        pb.xlower[ix_["SO"]] = 0.0
        pb.xlower[ix_["SS"]] = 0.0
        pb.xlower[ix_["IMPDER"]] = 100.0
        pb.xlower[ix_["IMPK"]] = 500.0
        pb.xlower[ix_["IMPFUS"]] = 500.0
        pb.xlower[ix_["QI"]] = 1000.0
        pb.xlower[ix_["PT"]] = 2.0
        pb.xlower[ix_["MV"]] = 2000.0
        pb.xlower[ix_["MC"]] = 3000.0
        pb.xlower[ix_["MD"]] = 5000.0
        pb.xlower[ix_["PD"]] = 0.2
        pb.xlower[ix_["NS"]] = 1.0
        pb.xlower[ix_["VS"]] = 0.0
        pb.xlower[ix_["CR"]] = 100.0
        pb.xlower[ix_["PM"]] = 4.0
        pb.xlower[ix_["DV"]] = 0.0
        pb.xlower[ix_["MZ"]] = 500.0
        pb.xlower[ix_["VN"]] = 10.0
        pb.xlower[ix_["QV"]] = 250.0
        pb.xlower[ix_["QF"]] = 750.0
        pb.xlower[ix_["IMPTRAIN"]] = 250.0
        pb.xlower[ix_["IMPMOT"]] = 10.0
        pb.xlower[ix_["IMPNMOT"]] = 35.0
        pb.xlower[ix_["IMPPET"]] = 100.0
        pb.xlower[ix_["IMPPIL"]] = 200.0
        pb.xlower[ix_["IMPCAN"]] = 120.0
        pb.xlower[ix_["IMPSNA"]] = 700.0
        pb.xlower[ix_["MS"]] = 100.0
        pb.xlower[ix_["EL"]] = 2.0
        pb.xlower[ix_["DE"]] = 0.0
        pb.xlower[ix_["DS"]] = 0.0
        pb.xlower[ix_["IMPVOIL"]] = 500.0
        pb.xupper[ix_["SR"]] = 150.0
        pb.xupper[ix_["LR"]] = 10.0
        pb.xupper[ix_["PK"]] = 10.0
        pb.xupper[ix_["EF"]] = 5.0
        pb.xupper[ix_["SX"]] = 120.0
        pb.xupper[ix_["LX"]] = 8.0
        pb.xupper[ix_["SD"]] = 20.0
        pb.xupper[ix_["SK"]] = 30.0
        pb.xupper[ix_["ST"]] = 500.0
        pb.xupper[ix_["SF"]] = 200.0
        pb.xupper[ix_["LF"]] = 20.0
        pb.xupper[ix_["AM"]] = 10.0
        pb.xupper[ix_["CA"]] = -0.001
        pb.xupper[ix_["CB"]] = 2.0
        pb.xupper[ix_["SO"]] = 1.0
        pb.xupper[ix_["SS"]] = 2.0
        pb.xupper[ix_["IMPDER"]] = 1000.0
        pb.xupper[ix_["IMPK"]] = 5000.0
        pb.xupper[ix_["IMPFUS"]] = 5000.0
        pb.xupper[ix_["QI"]] = 20000.0
        pb.xupper[ix_["PT"]] = 30.0
        pb.xupper[ix_["MV"]] = 20000.0
        pb.xupper[ix_["MC"]] = 30000.0
        pb.xupper[ix_["MD"]] = 50000.0
        pb.xupper[ix_["PD"]] = 0.8
        pb.xupper[ix_["NS"]] = 5.0
        pb.xupper[ix_["VS"]] = 20.0
        pb.xupper[ix_["CR"]] = 400.0
        pb.xupper[ix_["PM"]] = 15.0
        pb.xupper[ix_["DV"]] = 10.0
        pb.xupper[ix_["MZ"]] = 10000.0
        pb.xupper[ix_["VN"]] = 50.0
        pb.xupper[ix_["QV"]] = 5000.0
        pb.xupper[ix_["QF"]] = 15000.0
        pb.xupper[ix_["IMPTRAIN"]] = 3000.0
        pb.xupper[ix_["IMPMOT"]] = 5000.0
        pb.xupper[ix_["IMPNMOT"]] = 70.0
        pb.xupper[ix_["IMPPET"]] = 3000.0
        pb.xupper[ix_["IMPPIL"]] = 400.0
        pb.xupper[ix_["IMPCAN"]] = 240.0
        pb.xupper[ix_["IMPSNA"]] = 1900.0
        pb.xupper[ix_["MS"]] = 1000.0
        pb.xupper[ix_["EL"]] = 20.0
        pb.xupper[ix_["DE"]] = 1.0
        pb.xupper[ix_["DS"]] = 2.0
        pb.xupper[ix_["IMPVOIL"]] = 5000.0
        pb.xupper[ix_["NM"]] = 2.0
        pb.xlower[ix_["NM"]] = 1.0
        pb.xupper[ix_["NP"]] = 2.0
        pb.xlower[ix_["NP"]] = 1.0
        pb.xupper[ix_["NG"]] = 2.0
        pb.xlower[ix_["NG"]] = 1.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        pb.x0[ix_["SR"]] = Float64(2.7452e+01)
        pb.x0[ix_["LR"]] = Float64(1.5000)
        pb.x0[ix_["PK"]] = Float64(1.0000e+01)
        pb.x0[ix_["EF"]] = Float64(0.0000)
        pb.x0[ix_["SX"]] = Float64(1.9217e+01)
        pb.x0[ix_["LX"]] = Float64(1.5000)
        pb.x0[ix_["SD"]] = Float64(3.5688)
        pb.x0[ix_["SK"]] = Float64(4.0696)
        pb.x0[ix_["ST"]] = Float64(3.4315e+01)
        pb.x0[ix_["SF"]] = Float64(8.8025e+01)
        pb.x0[ix_["LF"]] = Float64(5.1306)
        pb.x0[ix_["AM"]] = Float64(0.0000)
        pb.x0[ix_["CA"]] = Float64(-1.4809e-01)
        pb.x0[ix_["CB"]] = Float64(7.5980e-01)
        pb.x0[ix_["SO"]] = Float64(0.0000)
        pb.x0[ix_["SS"]] = Float64(0.0000)
        pb.x0[ix_["IMPDER"]] = Float64(1.1470e+02)
        pb.x0[ix_["IMPK"]] = Float64(5.0000e+02)
        pb.x0[ix_["IMPFUS"]] = Float64(1.7605e+03)
        pb.x0[ix_["QI"]] = Float64(2.3256e+03)
        pb.x0[ix_["PT"]] = Float64(5.6788)
        pb.x0[ix_["MV"]] = Float64(1.4197e+04)
        pb.x0[ix_["MC"]] = Float64(1.2589e+04)
        pb.x0[ix_["MD"]] = Float64(2.8394e+04)
        pb.x0[ix_["PD"]] = Float64(2.0000e-01)
        pb.x0[ix_["NS"]] = Float64(1.0000)
        pb.x0[ix_["VS"]] = Float64(0.0000)
        pb.x0[ix_["CR"]] = Float64(1.0000e+02)
        pb.x0[ix_["PM"]] = Float64(1.5000e+01)
        pb.x0[ix_["DV"]] = Float64(0.0000)
        pb.x0[ix_["MZ"]] = Float64(5.0000e+02)
        pb.x0[ix_["VN"]] = Float64(1.0000e+01)
        pb.x0[ix_["QV"]] = Float64(8.1490e+02)
        pb.x0[ix_["QF"]] = Float64(3.1405e+03)
        pb.x0[ix_["IMPTRAIN"]] = Float64(1.9450e+03)
        pb.x0[ix_["IMPMOT"]] = Float64(1.9085e+02)
        pb.x0[ix_["IMPNMOT"]] = Float64(3.5000e+01)
        pb.x0[ix_["IMPPET"]] = Float64(1.0000e+02)
        pb.x0[ix_["IMPPIL"]] = Float64(2.0000e+02)
        pb.x0[ix_["IMPCAN"]] = Float64(1.2000e+02)
        pb.x0[ix_["IMPSNA"]] = Float64(7.0000e+02)
        pb.x0[ix_["MS"]] = Float64(1.0000e+03)
        pb.x0[ix_["EL"]] = Float64(4.9367)
        pb.x0[ix_["DE"]] = Float64(0.0000)
        pb.x0[ix_["DS"]] = Float64(0.0000)
        pb.x0[ix_["IMPVOIL"]] = Float64(5.0000e+03)
        pb.x0[ix_["NM"]] = Float64(1.0000)
        pb.x0[ix_["NP"]] = Float64(1.0000)
        pb.x0[ix_["NG"]] = Float64(1.0000)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "en2PR", iet_)
        loaset(elftv,it,1,"Y")
        loaset(elftv,it,2,"Z")
        it,iet_,_ = s2mpj_ii( "eQDdSQ", iet_)
        loaset(elftv,it,1,"W")
        loaset(elftv,it,2,"X")
        loaset(elftv,it,3,"Y")
        loaset(elftv,it,4,"Z")
        it,iet_,_ = s2mpj_ii( "en12", iet_)
        loaset(elftv,it,1,"Y")
        loaset(elftv,it,2,"Z")
        it,iet_,_ = s2mpj_ii( "en12d1", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        loaset(elftv,it,3,"Z")
        it,iet_,_ = s2mpj_ii( "eSQ", iet_)
        loaset(elftv,it,1,"Z")
        it,iet_,_ = s2mpj_ii( "eQT", iet_)
        loaset(elftv,it,1,"Y")
        loaset(elftv,it,2,"Z")
        it,iet_,_ = s2mpj_ii( "en1dLIN", iet_)
        loaset(elftv,it,1,"Y")
        loaset(elftv,it,2,"Z")
        it,iet_,_ = s2mpj_ii( "eSQRT", iet_)
        loaset(elftv,it,1,"Z")
        it,iet_,_ = s2mpj_ii( "eSURD", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        loaset(elftv,it,3,"Z")
        it,iet_,_ = s2mpj_ii( "eSQPRD", iet_)
        loaset(elftv,it,1,"Y")
        loaset(elftv,it,2,"Z")
        it,iet_,_ = s2mpj_ii( "eCBdSQQD", iet_)
        loaset(elftv,it,1,"W")
        loaset(elftv,it,2,"X")
        loaset(elftv,it,3,"Y")
        loaset(elftv,it,4,"Z")
        it,iet_,_ = s2mpj_ii( "eSREL", iet_)
        loaset(elftv,it,1,"V")
        loaset(elftv,it,2,"W")
        loaset(elftv,it,3,"X")
        loaset(elftv,it,4,"Y")
        loaset(elftv,it,5,"Z")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "EL1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en2PR")
        arrset(ielftype,ie,iet_["en2PR"])
        vname = "PK"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "SR"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EL2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eQDdSQ")
        arrset(ielftype,ie,iet_["eQDdSQ"])
        vname = "SS"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="W",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "SO"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "CB"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "LF"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EL3"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en12")
        arrset(ielftype,ie,iet_["en12"])
        vname = "EF"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "LF"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EL4"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en12d1")
        arrset(ielftype,ie,iet_["en12d1"])
        vname = "SO"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "CB"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "CA"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EL5"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "SD"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EL6"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "SK"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EL7"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "MV"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EL8"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en2PR")
        arrset(ielftype,ie,iet_["en2PR"])
        vname = "MD"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PD"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EL9"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eQT")
        arrset(ielftype,ie,iet_["eQT"])
        vname = "MZ"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "CR"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EL10"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en2PR")
        arrset(ielftype,ie,iet_["en2PR"])
        vname = "DV"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PT"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EL11"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en1dLIN")
        arrset(ielftype,ie,iet_["en1dLIN"])
        vname = "PT"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PM"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EL12"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQRT")
        arrset(ielftype,ie,iet_["eSQRT"])
        vname = "PT"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EL13"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en2PR")
        arrset(ielftype,ie,iet_["en2PR"])
        vname = "SR"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "NM"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EL14"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eQT")
        arrset(ielftype,ie,iet_["eQT"])
        vname = "MD"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "MS"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EL15"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSURD")
        arrset(ielftype,ie,iet_["eSURD"])
        vname = "SX"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "EL"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "LX"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EL16"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQPRD")
        arrset(ielftype,ie,iet_["eSQPRD"])
        vname = "DE"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PT"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EL17"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQPRD")
        arrset(ielftype,ie,iet_["eSQPRD"])
        vname = "DS"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PT"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EL18"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eCBdSQQD")
        arrset(ielftype,ie,iet_["eCBdSQQD"])
        vname = "VN"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="W",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "CA"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "LF"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "SO"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EL19"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSREL")
        arrset(ielftype,ie,iet_["eSREL"])
        vname = "SX"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "MC"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="W",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "LX"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "SR"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "EL"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gSQUARE",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["E4"]
        arrset(pbm.grftype,ig,"gSQUARE")
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EL1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-0.01))
        ig = ig_["E6"]
        arrset(pbm.grftype,ig,"gSQUARE")
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EL2"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        ig = ig_["E7"]
        arrset(pbm.grftype,ig,"gSQUARE")
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EL3"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(0.01))
        ig = ig_["E8"]
        arrset(pbm.grftype,ig,"gSQUARE")
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EL4"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-0.25))
        ig = ig_["E9"]
        arrset(pbm.grftype,ig,"gSQUARE")
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EL5"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-1.3))
        ig = ig_["E10"]
        arrset(pbm.grftype,ig,"gSQUARE")
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EL6"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(8.6))
        ig = ig_["E13"]
        arrset(pbm.grftype,ig,"gSQUARE")
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EL7"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["1/24000"]))
        ig = ig_["E14"]
        arrset(pbm.grftype,ig,"gSQUARE")
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EL8"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        ig = ig_["E16"]
        arrset(pbm.grftype,ig,"gSQUARE")
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EL9"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["EL10"])
        loaset(pbm.grelw,ig,posel,Float64(1.0))
        ig = ig_["E18"]
        arrset(pbm.grftype,ig,"gSQUARE")
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EL11"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-1000.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["EL12"])
        loaset(pbm.grelw,ig,posel,Float64(-12.0))
        ig = ig_["E26"]
        arrset(pbm.grftype,ig,"gSQUARE")
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EL13"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-1.25))
        ig = ig_["E27"]
        arrset(pbm.grftype,ig,"gSQUARE")
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EL14"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        ig = ig_["E28"]
        arrset(pbm.grftype,ig,"gSQUARE")
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EL15"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-2.4))
        ig = ig_["E29"]
        arrset(pbm.grftype,ig,"gSQUARE")
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EL16"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-0.785))
        ig = ig_["E30"]
        arrset(pbm.grftype,ig,"gSQUARE")
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EL17"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-0.785))
        ig = ig_["E31"]
        arrset(pbm.grftype,ig,"gSQUARE")
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EL18"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-2.0))
        ig = ig_["E32"]
        arrset(pbm.grftype,ig,"gSQUARE")
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EL19"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-1.15))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               9.46801297093018D+07
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
        pb.pbclass = "C-OLR2-RN-49-15"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm

# ***********************
#  SET UP THE FUNCTIONS *
# ***********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

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
                H_[1,2] = 1.0e0
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

    elseif action == "eQDdSQ"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        QD = EV_[1]-EV_[2]-EV_[3]*EV_[4]
        SQ = EV_[4]^2
        RSQ = 1.0e0/SQ
        QDOSQ = QD/SQ
        f_   = QDOSQ
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = RSQ
            g_[2] = -RSQ
            g_[3] = -1.0e0/EV_[4]
            g_[4] = -EV_[3]/SQ-2.0e0*QDOSQ/EV_[4]
            if nargout>2
                H_ = zeros(Float64,4,4)
                H_[1,4] = -2.0e0/(SQ*EV_[4])
                H_[4,1] = H_[1,4]
                H_[2,4] = 2.0e0/(SQ*EV_[4])
                H_[4,2] = H_[2,4]
                H_[3,4] = RSQ
                H_[4,3] = H_[3,4]
                H_[4,4] = (4.0e0*EV_[3])/(SQ*EV_[4])+6.0e0*QDOSQ/SQ
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "en12"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]/EV_[2]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 1.0/EV_[2]
            g_[2] = -EV_[1]/EV_[2]^2
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = -1.0e0/EV_[2]^2
                H_[2,1] = H_[1,2]
                H_[2,2] = (2.0e0*EV_[1])/EV_[2]^3
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "en12d1"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        ZSQ = EV_[3]^2
        YSQ = EV_[2]^2
        f_   = (EV_[1]*YSQ)/EV_[3]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = YSQ/EV_[3]
            g_[2] = (2.0e0*EV_[1]*EV_[2])/EV_[3]
            g_[3] = -(EV_[1]*YSQ)/ZSQ
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,2] = (2.0e0*EV_[2])/EV_[3]
                H_[2,1] = H_[1,2]
                H_[1,3] = -YSQ/ZSQ
                H_[3,1] = H_[1,3]
                H_[2,2] = (2.0e0*EV_[1])/EV_[3]
                H_[2,3] = -(2.0e0*EV_[1]*EV_[2])/ZSQ
                H_[3,2] = H_[2,3]
                H_[3,3] = (2.0e0*EV_[1]*YSQ)/(ZSQ*EV_[3])
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eSQ"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 2.0e0*EV_[1]
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 2.0e0
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eQT"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]/EV_[2]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 1.0/EV_[2]
            g_[2] = -EV_[1]/EV_[2]^2
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = -1.0e0/EV_[2]^2
                H_[2,1] = H_[1,2]
                H_[2,2] = (2.0e0*EV_[1])/EV_[2]^3
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "en1dLIN"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        LIN = EV_[2]+20.0e0
        f_   = EV_[1]/LIN
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 1.0/LIN
            g_[2] = -EV_[1]/LIN^2
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = -1.0e0/LIN^2
                H_[2,1] = H_[1,2]
                H_[2,2] = (2.0e0*EV_[1])/LIN^3
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eSQRT"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        RTZ = sqrt(EV_[1])
        f_   = RTZ
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 0.5e0/RTZ
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = -0.25e0/(EV_[1]*RTZ)
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eSURD"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        RTX = sqrt(EV_[1])
        RTZ = sqrt(EV_[3])
        XRTX = EV_[1]*sqrt(EV_[1])
        ZRTZ = EV_[3]*sqrt(EV_[3])
        f_   = XRTX*EV_[2]/RTZ
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 1.5e0*RTX*EV_[2]/RTZ
            g_[2] = XRTX/RTZ
            g_[3] = -(0.5e0*XRTX*EV_[2])/ZRTZ
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,1] = 0.75e0*EV_[2]/(RTZ*RTX)
                H_[1,2] = 1.5e0*RTX/RTZ
                H_[2,1] = H_[1,2]
                H_[1,3] = -(0.75e0*RTX*EV_[2])/ZRTZ
                H_[3,1] = H_[1,3]
                H_[2,3] = -(0.5e0*XRTX)/ZRTZ
                H_[3,2] = H_[2,3]
                H_[3,3] = (0.75e0*XRTX*EV_[2])/(ZRTZ*EV_[3])
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eSQPRD"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]^2*EV_[2]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 2.0e0*EV_[1]*EV_[2]
            g_[2] = EV_[1]^2
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = 2.0e0*EV_[2]
                H_[1,2] = 2.0e0*EV_[1]
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

    elseif action == "eCBdSQQD"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        YCB = EV_[3]^3
        CB = EV_[1]-EV_[2]*YCB
        TMZY = 3.0e0-EV_[4]*EV_[3]
        SQQD = EV_[3]^2*TMZY
        DYSQQD = 3.0e0*EV_[3]*(2.0e0-EV_[4]*EV_[3])
        f_   = CB/SQQD
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 1.0e0/SQQD
            g_[2] = -YCB/SQQD
            g_[3] = -3.0e0*EV_[2]*EV_[3]^2/SQQD-(CB*DYSQQD)/(SQQD*SQQD)
            g_[4] = (YCB*CB)/SQQD^2
            if nargout>2
                H_ = zeros(Float64,4,4)
                H_[1,3] = -DYSQQD/SQQD^2
                H_[3,1] = H_[1,3]
                H_[1,4] = YCB/SQQD^2
                H_[4,1] = H_[1,4]
                H_[2,3] = -3.0e0*(1.0e0-(TMZY-1.0e0)/TMZY)/TMZY
                H_[3,2] = H_[2,3]
                H_[2,4] = -(YCB^2*CB)/SQQD^2
                H_[4,2] = H_[2,4]
                H_[3,3] = (-6.0e0*EV_[2]*EV_[3]/SQQD+(3.0e0*EV_[2]*EV_[3]^2*DYSQQD)/
                     SQQD^2+(2.0e0*CB*DYSQQD^2)/SQQD^3-(6.0e0*CB*(1.0-EV_[4]*EV_[3]))/SQQD^2)
                H_[3,4] = (-(3.0e0*EV_[2]*EV_[3])/TMZY^2-(6.0e0*EV_[3]^4*(TMZY-1.0e0)*CB)/
                     SQQD^3+3.0e0*CB*EV_[3]^2/SQQD^2)
                H_[4,3] = H_[3,4]
                H_[4,4] = 2.0e0*(YCB^2*CB)/SQQD^3
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eSREL"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        SRLIN = 15.0e0+0.15e0*EV_[1]
        SRPD = EV_[2]*EV_[3]
        SRQD = 50.0e0*EV_[4]*EV_[5]
        SRQT = SRPD/SRQD
        SRRT = sqrt(SRQT)
        SRSUM = 15.0e0+0.3e0*EV_[1]
        SRVLN = EV_[1]*SRLIN
        f_   = EV_[1]*SRLIN*(SRQT*SRRT+8.0e0)
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = SRSUM*(SRQT*SRRT+8.0e0)
            g_[2] = 1.5e0*(SRVLN*SRRT*EV_[3]/SRQD)
            g_[3] = 1.5e0*(SRVLN*SRRT*EV_[2]/SRQD)
            g_[4] = -1.5e0*SRVLN*SRRT*SRQT/EV_[4]
            g_[5] = -1.5e0*SRVLN*SRRT*SRQT/EV_[5]
            if nargout>2
                H_ = zeros(Float64,5,5)
                H_[1,1] = 0.3e0*(SRQT*SRRT+8.0)
                H_[1,2] = 1.5e0*(SRSUM*SRRT*EV_[3]/SRQD)
                H_[2,1] = H_[1,2]
                H_[1,3] = 1.5e0*(SRSUM*SRRT*EV_[2]/SRQD)
                H_[3,1] = H_[1,3]
                H_[1,4] = -1.5e0*SRSUM*SRRT*SRQT/EV_[4]
                H_[4,1] = H_[1,4]
                H_[1,5] = -1.5e0*SRSUM*SRRT*SRQT/EV_[5]
                H_[5,1] = H_[1,5]
                H_[2,2] = (0.75e0*SRVLN*EV_[3]^2)/(SRQD^2*SRRT)
                H_[2,3] = SRVLN*((0.75e0*SRPD)/(SRQD^2*SRRT)+(1.5e0*SRRT)/SRQD)
                H_[3,2] = H_[2,3]
                H_[2,4] = -(SRVLN*2.25e0*SRRT*SRQT)/(EV_[2]*EV_[4])
                H_[4,2] = H_[2,4]
                H_[2,5] = -(SRVLN*2.25e0*SRRT*SRQT)/(EV_[2]*EV_[5])
                H_[5,2] = H_[2,5]
                H_[3,3] = (SRVLN*0.75e0*EV_[2]*EV_[2])/(SRRT*SRQD^2)
                H_[3,4] = -(SRVLN*2.25e0*SRRT*SRQT)/(EV_[3]*EV_[4])
                H_[4,3] = H_[3,4]
                H_[3,5] = -(SRVLN*2.25e0*SRRT*SRQT)/(EV_[3]*EV_[5])
                H_[5,3] = H_[3,5]
                H_[4,4] = (SRVLN*3.75e0*SRRT*SRQT)/EV_[4]^2
                H_[4,5] = (SRVLN*2.25e0*SRRT*SRQT)/(EV_[4]*EV_[5])
                H_[5,4] = H_[4,5]
                H_[5,5] = (SRVLN*3.75e0*SRRT*SRQT)/EV_[5]^2
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

    elseif action == "gSQUARE"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_= GVAR_*GVAR_
        if nargout>1
            g_ = GVAR_+GVAR_
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = 2.0e0
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

