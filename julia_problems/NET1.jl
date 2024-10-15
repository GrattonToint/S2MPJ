function NET1(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : NET1
#    *********
# 
#    A gas network problem for the south-east of England.
# 
#     SIF input: Sybille Schachler, Oxford, August 1992.
#    classification = "C-OOI2-RN-48-57"
# 
#    ...Problem size parameters
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "NET1"

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
        v_["NNOD"] = 22
        v_["NPIP"] = 17
        v_["NCMP"] = 3
        v_["NSRC"] = 2
        v_["CSTART"] = 18
        v_["CEND"] = 20
        v_["SSTART"] = 21
        v_["SEND"] = 22
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        iv,ix_,_ = s2mpj_ii("NOP17",ix_)
        arrset(pb.xnames,iv,"NOP17")
        iv,ix_,_ = s2mpj_ii("PFL11",ix_)
        arrset(pb.xnames,iv,"PFL11")
        iv,ix_,_ = s2mpj_ii("NOP9",ix_)
        arrset(pb.xnames,iv,"NOP9")
        iv,ix_,_ = s2mpj_ii("PFL10",ix_)
        arrset(pb.xnames,iv,"PFL10")
        iv,ix_,_ = s2mpj_ii("NOP16",ix_)
        arrset(pb.xnames,iv,"NOP16")
        iv,ix_,_ = s2mpj_ii("PFL16",ix_)
        arrset(pb.xnames,iv,"PFL16")
        iv,ix_,_ = s2mpj_ii("NOP19",ix_)
        arrset(pb.xnames,iv,"NOP19")
        iv,ix_,_ = s2mpj_ii("PFL17",ix_)
        arrset(pb.xnames,iv,"PFL17")
        iv,ix_,_ = s2mpj_ii("SFL22",ix_)
        arrset(pb.xnames,iv,"SFL22")
        iv,ix_,_ = s2mpj_ii("SBV22",ix_)
        arrset(pb.xnames,iv,"SBV22")
        iv,ix_,_ = s2mpj_ii("NOP18",ix_)
        arrset(pb.xnames,iv,"NOP18")
        iv,ix_,_ = s2mpj_ii("NOP4",ix_)
        arrset(pb.xnames,iv,"NOP4")
        iv,ix_,_ = s2mpj_ii("PFL13",ix_)
        arrset(pb.xnames,iv,"PFL13")
        iv,ix_,_ = s2mpj_ii("PFL5",ix_)
        arrset(pb.xnames,iv,"PFL5")
        iv,ix_,_ = s2mpj_ii("NOP11",ix_)
        arrset(pb.xnames,iv,"NOP11")
        iv,ix_,_ = s2mpj_ii("PFL8",ix_)
        arrset(pb.xnames,iv,"PFL8")
        iv,ix_,_ = s2mpj_ii("NOP6",ix_)
        arrset(pb.xnames,iv,"NOP6")
        iv,ix_,_ = s2mpj_ii("NOP12",ix_)
        arrset(pb.xnames,iv,"NOP12")
        iv,ix_,_ = s2mpj_ii("CFL19",ix_)
        arrset(pb.xnames,iv,"CFL19")
        iv,ix_,_ = s2mpj_ii("CBV19",ix_)
        arrset(pb.xnames,iv,"CBV19")
        iv,ix_,_ = s2mpj_ii("PFL7",ix_)
        arrset(pb.xnames,iv,"PFL7")
        iv,ix_,_ = s2mpj_ii("NOP5",ix_)
        arrset(pb.xnames,iv,"NOP5")
        iv,ix_,_ = s2mpj_ii("PFL6",ix_)
        arrset(pb.xnames,iv,"PFL6")
        iv,ix_,_ = s2mpj_ii("NOP8",ix_)
        arrset(pb.xnames,iv,"NOP8")
        iv,ix_,_ = s2mpj_ii("CFL20",ix_)
        arrset(pb.xnames,iv,"CFL20")
        iv,ix_,_ = s2mpj_ii("CBV20",ix_)
        arrset(pb.xnames,iv,"CBV20")
        iv,ix_,_ = s2mpj_ii("NOP7",ix_)
        arrset(pb.xnames,iv,"NOP7")
        iv,ix_,_ = s2mpj_ii("PFL9",ix_)
        arrset(pb.xnames,iv,"PFL9")
        iv,ix_,_ = s2mpj_ii("NOP21",ix_)
        arrset(pb.xnames,iv,"NOP21")
        iv,ix_,_ = s2mpj_ii("PFL2",ix_)
        arrset(pb.xnames,iv,"PFL2")
        iv,ix_,_ = s2mpj_ii("SFL21",ix_)
        arrset(pb.xnames,iv,"SFL21")
        iv,ix_,_ = s2mpj_ii("SBV21",ix_)
        arrset(pb.xnames,iv,"SBV21")
        iv,ix_,_ = s2mpj_ii("NOP1",ix_)
        arrset(pb.xnames,iv,"NOP1")
        iv,ix_,_ = s2mpj_ii("PFL1",ix_)
        arrset(pb.xnames,iv,"PFL1")
        iv,ix_,_ = s2mpj_ii("NOP14",ix_)
        arrset(pb.xnames,iv,"NOP14")
        iv,ix_,_ = s2mpj_ii("PFL12",ix_)
        arrset(pb.xnames,iv,"PFL12")
        iv,ix_,_ = s2mpj_ii("NOP10",ix_)
        arrset(pb.xnames,iv,"NOP10")
        iv,ix_,_ = s2mpj_ii("PFL3",ix_)
        arrset(pb.xnames,iv,"PFL3")
        iv,ix_,_ = s2mpj_ii("NOP2",ix_)
        arrset(pb.xnames,iv,"NOP2")
        iv,ix_,_ = s2mpj_ii("CFL18",ix_)
        arrset(pb.xnames,iv,"CFL18")
        iv,ix_,_ = s2mpj_ii("CBV18",ix_)
        arrset(pb.xnames,iv,"CBV18")
        iv,ix_,_ = s2mpj_ii("NOP3",ix_)
        arrset(pb.xnames,iv,"NOP3")
        iv,ix_,_ = s2mpj_ii("PFL4",ix_)
        arrset(pb.xnames,iv,"PFL4")
        iv,ix_,_ = s2mpj_ii("NOP15",ix_)
        arrset(pb.xnames,iv,"NOP15")
        iv,ix_,_ = s2mpj_ii("PFL15",ix_)
        arrset(pb.xnames,iv,"PFL15")
        iv,ix_,_ = s2mpj_ii("NOP20",ix_)
        arrset(pb.xnames,iv,"NOP20")
        iv,ix_,_ = s2mpj_ii("PFL14",ix_)
        arrset(pb.xnames,iv,"PFL14")
        iv,ix_,_ = s2mpj_ii("NOP13",ix_)
        arrset(pb.xnames,iv,"NOP13")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("MBE1",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"MBE1")
        iv = ix_["PFL1"]
        pbm.A[ig,iv] += Float64(-1)
        iv = ix_["PFL2"]
        pbm.A[ig,iv] += Float64(-1)
        iv = ix_["SFL21"]
        pbm.A[ig,iv] += Float64(1)
        ig,ig_,_ = s2mpj_ii("MBE2",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"MBE2")
        iv = ix_["PFL3"]
        pbm.A[ig,iv] += Float64(-1)
        iv = ix_["CFL18"]
        pbm.A[ig,iv] += Float64(-1)
        ig,ig_,_ = s2mpj_ii("MBE3",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"MBE3")
        iv = ix_["PFL4"]
        pbm.A[ig,iv] += Float64(-1)
        iv = ix_["CFL18"]
        pbm.A[ig,iv] += Float64(1)
        ig,ig_,_ = s2mpj_ii("MBE4",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"MBE4")
        iv = ix_["PFL5"]
        pbm.A[ig,iv] += Float64(-1)
        iv = ix_["SFL22"]
        pbm.A[ig,iv] += Float64(1)
        ig,ig_,_ = s2mpj_ii("MBE5",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"MBE5")
        iv = ix_["PFL6"]
        pbm.A[ig,iv] += Float64(-1)
        iv = ix_["PFL7"]
        pbm.A[ig,iv] += Float64(-1)
        iv = ix_["CFL19"]
        pbm.A[ig,iv] += Float64(-1)
        ig,ig_,_ = s2mpj_ii("MBE6",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"MBE6")
        iv = ix_["PFL8"]
        pbm.A[ig,iv] += Float64(-1)
        iv = ix_["CFL19"]
        pbm.A[ig,iv] += Float64(1)
        ig,ig_,_ = s2mpj_ii("MBE7",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"MBE7")
        iv = ix_["PFL9"]
        pbm.A[ig,iv] += Float64(-1)
        iv = ix_["CFL20"]
        pbm.A[ig,iv] += Float64(-1)
        ig,ig_,_ = s2mpj_ii("MBE8",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"MBE8")
        iv = ix_["PFL6"]
        pbm.A[ig,iv] += Float64(1)
        iv = ix_["CFL20"]
        pbm.A[ig,iv] += Float64(1)
        ig,ig_,_ = s2mpj_ii("MBE9",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"MBE9")
        iv = ix_["PFL10"]
        pbm.A[ig,iv] += Float64(-1)
        iv = ix_["PFL11"]
        pbm.A[ig,iv] += Float64(-1)
        ig,ig_,_ = s2mpj_ii("MBE10",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"MBE10")
        iv = ix_["PFL3"]
        pbm.A[ig,iv] += Float64(1)
        iv = ix_["PFL12"]
        pbm.A[ig,iv] += Float64(-1)
        ig,ig_,_ = s2mpj_ii("MBE11",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"MBE11")
        iv = ix_["PFL5"]
        pbm.A[ig,iv] += Float64(1)
        iv = ix_["PFL8"]
        pbm.A[ig,iv] += Float64(1)
        iv = ix_["PFL13"]
        pbm.A[ig,iv] += Float64(-1)
        ig,ig_,_ = s2mpj_ii("MBE12",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"MBE12")
        iv = ix_["PFL7"]
        pbm.A[ig,iv] += Float64(1)
        ig,ig_,_ = s2mpj_ii("MBE13",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"MBE13")
        iv = ix_["PFL14"]
        pbm.A[ig,iv] += Float64(-1)
        ig,ig_,_ = s2mpj_ii("MBE14",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"MBE14")
        iv = ix_["PFL1"]
        pbm.A[ig,iv] += Float64(1)
        iv = ix_["PFL12"]
        pbm.A[ig,iv] += Float64(1)
        ig,ig_,_ = s2mpj_ii("MBE15",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"MBE15")
        iv = ix_["PFL4"]
        pbm.A[ig,iv] += Float64(1)
        iv = ix_["PFL15"]
        pbm.A[ig,iv] += Float64(-1)
        ig,ig_,_ = s2mpj_ii("MBE16",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"MBE16")
        iv = ix_["PFL10"]
        pbm.A[ig,iv] += Float64(1)
        iv = ix_["PFL16"]
        pbm.A[ig,iv] += Float64(-1)
        ig,ig_,_ = s2mpj_ii("MBE17",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"MBE17")
        iv = ix_["PFL11"]
        pbm.A[ig,iv] += Float64(1)
        ig,ig_,_ = s2mpj_ii("MBE18",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"MBE18")
        iv = ix_["PFL13"]
        pbm.A[ig,iv] += Float64(1)
        iv = ix_["PFL17"]
        pbm.A[ig,iv] += Float64(-1)
        ig,ig_,_ = s2mpj_ii("MBE19",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"MBE19")
        iv = ix_["PFL16"]
        pbm.A[ig,iv] += Float64(1)
        iv = ix_["PFL17"]
        pbm.A[ig,iv] += Float64(1)
        ig,ig_,_ = s2mpj_ii("MBE20",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"MBE20")
        iv = ix_["PFL14"]
        pbm.A[ig,iv] += Float64(1)
        iv = ix_["PFL15"]
        pbm.A[ig,iv] += Float64(1)
        ig,ig_,_ = s2mpj_ii("MBE21",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"MBE21")
        iv = ix_["PFL2"]
        pbm.A[ig,iv] += Float64(1)
        iv = ix_["PFL9"]
        pbm.A[ig,iv] += Float64(1)
        ig,ig_,_ = s2mpj_ii("MCR18",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MCR18")
        iv = ix_["NOP3"]
        pbm.A[ig,iv] += Float64(1.00000000)
        iv = ix_["NOP2"]
        pbm.A[ig,iv] += Float64(-1.40000000)
        iv = ix_["CBV18"]
        pbm.A[ig,iv] += Float64(0.00000)
        ig,ig_,_ = s2mpj_ii("MCR19",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MCR19")
        iv = ix_["NOP6"]
        pbm.A[ig,iv] += Float64(1.00000000)
        iv = ix_["NOP5"]
        pbm.A[ig,iv] += Float64(-1.40000000)
        iv = ix_["CBV19"]
        pbm.A[ig,iv] += Float64(0.00000)
        ig,ig_,_ = s2mpj_ii("MCR20",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"MCR20")
        iv = ix_["NOP8"]
        pbm.A[ig,iv] += Float64(1.00000000)
        iv = ix_["NOP7"]
        pbm.A[ig,iv] += Float64(-1.40000000)
        iv = ix_["CBV20"]
        pbm.A[ig,iv] += Float64(0.00000)
        ig,ig_,_ = s2mpj_ii("CLF18",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CLF18")
        iv = ix_["CBV18"]
        pbm.A[ig,iv] += Float64(1.00000e+04)
        iv = ix_["CFL18"]
        pbm.A[ig,iv] += Float64(-1.00000000)
        ig,ig_,_ = s2mpj_ii("CLF19",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CLF19")
        iv = ix_["CBV19"]
        pbm.A[ig,iv] += Float64(1.00000e+04)
        iv = ix_["CFL19"]
        pbm.A[ig,iv] += Float64(-1.00000000)
        ig,ig_,_ = s2mpj_ii("CLF20",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CLF20")
        iv = ix_["CBV20"]
        pbm.A[ig,iv] += Float64(1.00000e+04)
        iv = ix_["CFL20"]
        pbm.A[ig,iv] += Float64(-1.00000000)
        ig,ig_,_ = s2mpj_ii("SLF21",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"SLF21")
        iv = ix_["SBV21"]
        pbm.A[ig,iv] += Float64(0.00000)
        iv = ix_["SFL21"]
        pbm.A[ig,iv] += Float64(-1.00000000)
        ig,ig_,_ = s2mpj_ii("SUF21",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"SUF21")
        iv = ix_["SBV21"]
        pbm.A[ig,iv] += Float64(-3.00000e+03)
        iv = ix_["SFL21"]
        pbm.A[ig,iv] += Float64(1.00000000)
        ig,ig_,_ = s2mpj_ii("SLF22",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"SLF22")
        iv = ix_["SBV22"]
        pbm.A[ig,iv] += Float64(0.00000)
        iv = ix_["SFL22"]
        pbm.A[ig,iv] += Float64(-1.00000000)
        ig,ig_,_ = s2mpj_ii("SUF22",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"SUF22")
        iv = ix_["SBV22"]
        pbm.A[ig,iv] += Float64(-1.06000e+02)
        iv = ix_["SFL22"]
        pbm.A[ig,iv] += Float64(1.00000000)
        ig,ig_,_ = s2mpj_ii("CLP18",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CLP18")
        iv = ix_["NOP3"]
        pbm.A[ig,iv] += Float64(1)
        iv = ix_["NOP2"]
        pbm.A[ig,iv] += Float64(-1)
        iv = ix_["CBV18"]
        pbm.A[ig,iv] += Float64(-4.72000e+02)
        ig,ig_,_ = s2mpj_ii("CUP18",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CUP18")
        iv = ix_["NOP3"]
        pbm.A[ig,iv] += Float64(-1)
        iv = ix_["NOP2"]
        pbm.A[ig,iv] += Float64(1)
        ig,ig_,_ = s2mpj_ii("CLP19",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CLP19")
        iv = ix_["NOP6"]
        pbm.A[ig,iv] += Float64(1)
        iv = ix_["NOP5"]
        pbm.A[ig,iv] += Float64(-1)
        iv = ix_["CBV19"]
        pbm.A[ig,iv] += Float64(-3.45000e+02)
        ig,ig_,_ = s2mpj_ii("CUP19",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CUP19")
        iv = ix_["NOP6"]
        pbm.A[ig,iv] += Float64(-1)
        iv = ix_["NOP5"]
        pbm.A[ig,iv] += Float64(1)
        ig,ig_,_ = s2mpj_ii("CLP20",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CLP20")
        iv = ix_["NOP8"]
        pbm.A[ig,iv] += Float64(1)
        iv = ix_["NOP7"]
        pbm.A[ig,iv] += Float64(-1)
        iv = ix_["CBV20"]
        pbm.A[ig,iv] += Float64(-5.75000e+02)
        ig,ig_,_ = s2mpj_ii("CUP20",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CUP20")
        iv = ix_["NOP8"]
        pbm.A[ig,iv] += Float64(-1)
        iv = ix_["NOP7"]
        pbm.A[ig,iv] += Float64(1)
        for i = Int64(v_["1"]):Int64(v_["NPIP"])
            ig,ig_,_ = s2mpj_ii("PDE"*string(i),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"PDE"*string(i))
            arrset(pbm.gscale,ig,Float64(1.00000e+03))
        end
        for i = Int64(v_["CSTART"]):Int64(v_["CEND"])
            ig,ig_,_ = s2mpj_ii("HPCON"*string(i),ig_)
            arrset(gtype,ig,"<=")
            arrset(pb.cnames,ig,"HPCON"*string(i))
            arrset(pbm.gscale,ig,Float64(70.00000000))
        end
        for i = Int64(v_["CSTART"]):Int64(v_["CEND"])
            ig,ig_,_ = s2mpj_ii("HPOBJ"*string(i),ig_)
            arrset(gtype,ig,"<>")
            arrset(pbm.gscale,ig,Float64(0.03500000))
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
        pbm.gconst[ig_["MBE1"]] = Float64(6.65680000)
        pbm.gconst[ig_["MBE4"]] = Float64(1.96100000)
        pbm.gconst[ig_["MBE9"]] = Float64(3.72060e+02)
        pbm.gconst[ig_["MBE10"]] = Float64(47.17000000)
        pbm.gconst[ig_["MBE11"]] = Float64(1.60060e+02)
        pbm.gconst[ig_["MBE12"]] = Float64(4.25060e+02)
        pbm.gconst[ig_["MBE13"]] = Float64(5.30000e+02)
        pbm.gconst[ig_["MBE14"]] = Float64(24.16800000)
        pbm.gconst[ig_["MBE15"]] = Float64(2.54400000)
        pbm.gconst[ig_["MBE16"]] = Float64(89.14600000)
        pbm.gconst[ig_["MBE17"]] = Float64(4.92900e+02)
        pbm.gconst[ig_["MBE20"]] = Float64(4.64280e+02)
        pbm.gconst[ig_["MBE21"]] = Float64(1.48400e+02)
        pbm.gconst[ig_["CLF18"]] = Float64(1.00000e+04)
        pbm.gconst[ig_["CLF19"]] = Float64(1.00000e+04)
        pbm.gconst[ig_["CLF20"]] = Float64(1.00000e+04)
        pbm.gconst[ig_["HPCON18"]] = Float64(2.07000e+04)
        pbm.gconst[ig_["HPCON19"]] = Float64(2.07000e+04)
        pbm.gconst[ig_["HPCON20"]] = Float64(4.14000e+04)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        pb.xlower[ix_["SFL21"]] = 0.00000
        pb.xupper[ix_["SFL21"]] = 3.00000e+03
        pb.xlower[ix_["SFL22"]] = 0.00000
        pb.xupper[ix_["SFL22"]] = 1.06000e+02
        pb.xlower[ix_["NOP1"]] = 5.00000e+02
        pb.xupper[ix_["NOP1"]] = 1.01500e+03
        pb.xlower[ix_["NOP2"]] = 5.00000e+02
        pb.xupper[ix_["NOP2"]] = 1.10000e+03
        pb.xlower[ix_["NOP3"]] = 5.00000e+02
        pb.xupper[ix_["NOP3"]] = 9.72000e+02
        pb.xlower[ix_["NOP4"]] = 5.00000e+02
        pb.xupper[ix_["NOP4"]] = 1.10000e+03
        pb.xlower[ix_["NOP5"]] = 5.00000e+02
        pb.xupper[ix_["NOP5"]] = 1.10000e+03
        pb.xlower[ix_["NOP6"]] = 5.00000e+02
        pb.xupper[ix_["NOP6"]] = 8.45000e+02
        pb.xlower[ix_["NOP7"]] = 5.00000e+02
        pb.xupper[ix_["NOP7"]] = 1.10000e+03
        pb.xlower[ix_["NOP8"]] = 5.00000e+02
        pb.xupper[ix_["NOP8"]] = 1.07500e+03
        pb.xlower[ix_["NOP9"]] = 5.00000e+02
        pb.xupper[ix_["NOP9"]] = 1.10000e+03
        pb.xlower[ix_["NOP10"]] = 5.00000e+02
        pb.xupper[ix_["NOP10"]] = 1.10000e+03
        pb.xlower[ix_["NOP11"]] = 5.00000e+02
        pb.xupper[ix_["NOP11"]] = 1.10000e+03
        pb.xlower[ix_["NOP12"]] = 5.00000e+02
        pb.xupper[ix_["NOP12"]] = 1.10000e+03
        pb.xlower[ix_["NOP13"]] = 5.80000e+02
        pb.xupper[ix_["NOP13"]] = 1.10000e+03
        pb.xlower[ix_["NOP14"]] = 5.00000e+02
        pb.xupper[ix_["NOP14"]] = 1.10000e+03
        pb.xlower[ix_["NOP15"]] = 5.00000e+02
        pb.xupper[ix_["NOP15"]] = 1.10000e+03
        pb.xlower[ix_["NOP16"]] = 5.00000e+02
        pb.xupper[ix_["NOP16"]] = 1.10000e+03
        pb.xlower[ix_["NOP17"]] = 5.00000e+02
        pb.xupper[ix_["NOP17"]] = 1.10000e+03
        pb.xlower[ix_["NOP18"]] = 5.00000e+02
        pb.xupper[ix_["NOP18"]] = 1.10000e+03
        pb.xlower[ix_["NOP19"]] = 5.00000e+02
        pb.xupper[ix_["NOP19"]] = 1.10000e+03
        pb.xlower[ix_["NOP20"]] = 5.00000e+02
        pb.xupper[ix_["NOP20"]] = 1.10000e+03
        pb.xlower[ix_["NOP21"]] = 5.00000e+02
        pb.xupper[ix_["NOP21"]] = 1.10000e+03
        pb.xlower[ix_["CBV18"]] = 0
        pb.xupper[ix_["CBV18"]] = 0
        pb.xlower[ix_["CBV19"]] = 1
        pb.xupper[ix_["CBV19"]] = 1
        pb.xlower[ix_["CBV20"]] = 1
        pb.xupper[ix_["CBV20"]] = 1
        pb.xlower[ix_["SBV21"]] = 1
        pb.xupper[ix_["SBV21"]] = 1
        pb.xlower[ix_["SBV22"]] = 1
        pb.xupper[ix_["SBV22"]] = 1
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(5.00000e+02),pb.n)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eA0PANHAN", iet_)
        loaset(elftv,it,1,"PIN")
        loaset(elftv,it,2,"POUT")
        loaset(elftv,it,3,"FLOW")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"PIPRES")
        it,iet_,_ = s2mpj_ii( "eA1MAXHP", iet_)
        loaset(elftv,it,1,"PIN")
        loaset(elftv,it,2,"POUT")
        loaset(elftv,it,3,"FLOW")
        loaset(elftv,it,4,"CBV")
        loaset(elftp,it,1,"IPL")
        loaset(elftp,it,2,"OPL")
        it,iet_,_ = s2mpj_ii( "eA2HPFUN", iet_)
        loaset(elftv,it,1,"PIN")
        loaset(elftv,it,2,"POUT")
        loaset(elftv,it,3,"FLOW")
        loaset(elftv,it,4,"CBV")
        loaset(elftp,it,1,"IPL")
        loaset(elftp,it,2,"OPL")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "PANH1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eA0PANHAN")
        arrset(ielftype,ie,iet_["eA0PANHAN"])
        vname = "NOP1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="PIN",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "NOP14"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="POUT",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PFL1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="FLOW",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="PIPRES",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.62131268))
        ename = "PANH2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eA0PANHAN")
        arrset(ielftype,ie,iet_["eA0PANHAN"])
        vname = "NOP1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="PIN",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "NOP21"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="POUT",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PFL2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="FLOW",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="PIPRES",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.31605264))
        ename = "PANH3"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eA0PANHAN")
        arrset(ielftype,ie,iet_["eA0PANHAN"])
        vname = "NOP2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="PIN",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "NOP10"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="POUT",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PFL3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="FLOW",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="PIPRES",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.13104611))
        ename = "PANH4"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eA0PANHAN")
        arrset(ielftype,ie,iet_["eA0PANHAN"])
        vname = "NOP3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="PIN",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "NOP15"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="POUT",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PFL4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="FLOW",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="PIPRES",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.12796251))
        ename = "PANH5"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eA0PANHAN")
        arrset(ielftype,ie,iet_["eA0PANHAN"])
        vname = "NOP4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="PIN",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "NOP11"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="POUT",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PFL5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="FLOW",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="PIPRES",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(3.78624623))
        ename = "PANH6"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eA0PANHAN")
        arrset(ielftype,ie,iet_["eA0PANHAN"])
        vname = "NOP5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="PIN",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "NOP8"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="POUT",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PFL6"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="FLOW",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="PIPRES",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.84948702))
        ename = "PANH7"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eA0PANHAN")
        arrset(ielftype,ie,iet_["eA0PANHAN"])
        vname = "NOP5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="PIN",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "NOP12"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="POUT",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PFL7"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="FLOW",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="PIPRES",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(2.13696026))
        ename = "PANH8"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eA0PANHAN")
        arrset(ielftype,ie,iet_["eA0PANHAN"])
        vname = "NOP6"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="PIN",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "NOP11"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="POUT",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PFL8"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="FLOW",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="PIPRES",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.25900862))
        ename = "PANH9"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eA0PANHAN")
        arrset(ielftype,ie,iet_["eA0PANHAN"])
        vname = "NOP7"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="PIN",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "NOP21"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="POUT",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PFL9"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="FLOW",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="PIPRES",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.32838618))
        ename = "PANH10"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eA0PANHAN")
        arrset(ielftype,ie,iet_["eA0PANHAN"])
        vname = "NOP9"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="PIN",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "NOP16"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="POUT",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PFL10"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="FLOW",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="PIPRES",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.33657520))
        ename = "PANH11"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eA0PANHAN")
        arrset(ielftype,ie,iet_["eA0PANHAN"])
        vname = "NOP9"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="PIN",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "NOP17"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="POUT",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PFL11"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="FLOW",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="PIPRES",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.61512113))
        ename = "PANH12"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eA0PANHAN")
        arrset(ielftype,ie,iet_["eA0PANHAN"])
        vname = "NOP10"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="PIN",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "NOP14"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="POUT",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PFL12"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="FLOW",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="PIPRES",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.51339271))
        ename = "PANH13"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eA0PANHAN")
        arrset(ielftype,ie,iet_["eA0PANHAN"])
        vname = "NOP11"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="PIN",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "NOP18"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="POUT",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PFL13"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="FLOW",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="PIPRES",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.20890923))
        ename = "PANH14"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eA0PANHAN")
        arrset(ielftype,ie,iet_["eA0PANHAN"])
        vname = "NOP13"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="PIN",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "NOP20"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="POUT",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PFL14"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="FLOW",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="PIPRES",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.15474706))
        ename = "PANH15"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eA0PANHAN")
        arrset(ielftype,ie,iet_["eA0PANHAN"])
        vname = "NOP15"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="PIN",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "NOP20"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="POUT",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PFL15"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="FLOW",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="PIPRES",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.26980036))
        ename = "PANH16"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eA0PANHAN")
        arrset(ielftype,ie,iet_["eA0PANHAN"])
        vname = "NOP16"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="PIN",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "NOP19"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="POUT",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PFL16"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="FLOW",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="PIPRES",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.04255562))
        ename = "PANH17"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eA0PANHAN")
        arrset(ielftype,ie,iet_["eA0PANHAN"])
        vname = "NOP18"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="PIN",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "NOP19"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="POUT",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PFL17"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="FLOW",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="PIPRES",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.12570329))
        ename = "HPMAX18"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eA1MAXHP")
        arrset(ielftype,ie,iet_["eA1MAXHP"])
        vname = "NOP2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="PIN",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "NOP3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="POUT",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "CFL18"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="FLOW",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "CBV18"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="CBV",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="IPL",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.00000))
        posep = findfirst(x->x=="OPL",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.00000))
        ename = "HPMAX19"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eA1MAXHP")
        arrset(ielftype,ie,iet_["eA1MAXHP"])
        vname = "NOP5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="PIN",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "NOP6"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="POUT",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "CFL19"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="FLOW",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "CBV19"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="CBV",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="IPL",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.00000))
        posep = findfirst(x->x=="OPL",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.00000))
        ename = "HPMAX20"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eA1MAXHP")
        arrset(ielftype,ie,iet_["eA1MAXHP"])
        vname = "NOP7"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="PIN",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "NOP8"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="POUT",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "CFL20"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="FLOW",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "CBV20"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="CBV",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="IPL",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.00000))
        posep = findfirst(x->x=="OPL",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.00000))
        ename = "HPFUN18"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eA2HPFUN")
        arrset(ielftype,ie,iet_["eA2HPFUN"])
        vname = "NOP2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="PIN",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "NOP3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="POUT",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "CFL18"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="FLOW",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "CBV18"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="CBV",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="IPL",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.00000))
        posep = findfirst(x->x=="OPL",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.00000))
        ename = "HPFUN19"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eA2HPFUN")
        arrset(ielftype,ie,iet_["eA2HPFUN"])
        vname = "NOP5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="PIN",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "NOP6"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="POUT",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "CFL19"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="FLOW",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "CBV19"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="CBV",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="IPL",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.00000))
        posep = findfirst(x->x=="OPL",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.00000))
        ename = "HPFUN20"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eA2HPFUN")
        arrset(ielftype,ie,iet_["eA2HPFUN"])
        vname = "NOP7"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="PIN",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "NOP8"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="POUT",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "CFL20"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="FLOW",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "CBV20"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(5.00000e+02))
        posev = findfirst(x->x=="CBV",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="IPL",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.00000))
        posep = findfirst(x->x=="OPL",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.00000))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for i = Int64(v_["1"]):Int64(v_["NPIP"])
            ig = ig_["PDE"*string(i)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["PANH"*string(i)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        for i = Int64(v_["CSTART"]):Int64(v_["CEND"])
            ig = ig_["HPCON"*string(i)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["HPMAX"*string(i)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        for i = Int64(v_["CSTART"]):Int64(v_["CEND"])
            ig = ig_["HPOBJ"*string(i)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["HPFUN"*string(i)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
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
        pb.pbclass = "C-OOI2-RN-48-57"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eA0PANHAN"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        A0FLEX = 1.8539e0
        A0FGT0 = EV_[3]>=0.0e0
        if A0FGT0
            A0HFLO = -pbm.elpar[iel_][1]*A0FLEX*(A0FLEX-1.0e0)*EV_[3]^(A0FLEX-2.0e0)
        end
        if !A0FGT0
            A0HFLO  = (
              A0FLEX*(A0FLEX-1.0e0)*pbm.elpar[iel_][1]*abs(EV_[3])^(A0FLEX-2.0e0))
        end
        f_    = (
              EV_[1]*EV_[1]-EV_[2]*EV_[2]-pbm.elpar[iel_][1]*EV_[3]*abs(EV_[3])^(A0FLEX-1.0e0))
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[2] = -2.0e0*EV_[2]
            g_[1] = 2.0e0*EV_[1]
            g_[3] = -pbm.elpar[iel_][1]*A0FLEX*abs(EV_[3])^(A0FLEX-1.0e0)
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[2,2] = -2.0e0
                H_[1,1] = 2.0e0
                H_[3,3] = A0HFLO
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eA1MAXHP"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        A1BETA = 0.23077e0
        A1HFAC = 203.712e0
        A1PSUC = EV_[1]-pbm.elpar[iel_][1]*EV_[4]
        A1PDIS = EV_[2]+pbm.elpar[iel_][2]*EV_[4]
        A1CRB = (A1PDIS/A1PSUC)^A1BETA
        A1PROD = A1BETA*A1HFAC*A1CRB*EV_[3]
        A1GPIN = -A1PROD/A1PSUC
        A1GPOU = A1PROD/A1PDIS
        A1GFLO = A1HFAC*(A1CRB-1.0e0)
        A1GCBV = -pbm.elpar[iel_][1]*A1GPIN+pbm.elpar[iel_][2]*A1GPOU
        A1HII = A1PROD*(A1BETA+1.0e0)/(A1PSUC^2)
        A1HIO = -A1PROD*A1BETA/(A1PSUC*A1PDIS)
        A1HOO = A1PROD*(A1BETA-1.0e0)/(A1PDIS^2)
        A1HIC = -pbm.elpar[iel_][1]*A1HII+pbm.elpar[iel_][2]*A1HIO
        A1HOC = -pbm.elpar[iel_][1]*A1HIO+pbm.elpar[iel_][2]*A1HOO
        f_   = EV_[3]*A1GFLO
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = A1GPIN
            g_[2] = A1GPOU
            g_[3] = A1GFLO
            g_[4] = A1GCBV
            if nargout>2
                H_ = zeros(Float64,4,4)
                H_[1,1] = A1HII
                H_[1,2] = A1HIO
                H_[2,1] = H_[1,2]
                H_[2,2] = A1HOO
                H_[1,3] = A1GPIN/EV_[3]
                H_[3,1] = H_[1,3]
                H_[2,3] = A1GPOU/EV_[3]
                H_[3,2] = H_[2,3]
                H_[1,4] = A1HIC
                H_[4,1] = H_[1,4]
                H_[2,4] = A1HOC
                H_[4,2] = H_[2,4]
                H_[3,4] = A1GCBV/EV_[3]
                H_[4,3] = H_[3,4]
                H_[4,4] = -pbm.elpar[iel_][1]*A1HIC+pbm.elpar[iel_][2]*A1HOC
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eA2HPFUN"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        A2BETA = 0.23077e0
        A2HFAC = 203.712e0
        A2PSUC = EV_[1]-pbm.elpar[iel_][1]*EV_[4]
        A2PDIS = EV_[2]+pbm.elpar[iel_][2]*EV_[4]
        A2CRB = (A2PDIS/A2PSUC)^A2BETA
        A2PROD = A2BETA*A2HFAC*A2CRB*EV_[3]
        A2GPIN = -A2PROD/A2PSUC
        A2GPOU = A2PROD/A2PDIS
        A2GFLO = A2HFAC*(A2CRB-1.0e0)
        A2GCBV = -pbm.elpar[iel_][1]*A2GPIN+pbm.elpar[iel_][2]*A2GPOU
        A2HII = A2PROD*(A2BETA+1.0e0)/(A2PSUC^2)
        A2HIO = -A2PROD*A2BETA/(A2PSUC*A2PDIS)
        A2HOO = A2PROD*(A2BETA-1.0e0)/(A2PDIS^2)
        A2HIC = -pbm.elpar[iel_][1]*A2HII+pbm.elpar[iel_][2]*A2HIO
        A2HOC = -pbm.elpar[iel_][1]*A2HIO+pbm.elpar[iel_][2]*A2HOO
        f_   = EV_[3]*A2GFLO
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = A2GPIN
            g_[2] = A2GPOU
            g_[3] = A2GFLO
            g_[4] = A2GCBV
            if nargout>2
                H_ = zeros(Float64,4,4)
                H_[1,1] = A2HII
                H_[1,2] = A2HIO
                H_[2,1] = H_[1,2]
                H_[2,2] = A2HOO
                H_[1,3] = A2GPIN/EV_[3]
                H_[3,1] = H_[1,3]
                H_[2,3] = A2GPOU/EV_[3]
                H_[3,2] = H_[2,3]
                H_[1,4] = A2HIC
                H_[4,1] = H_[1,4]
                H_[2,4] = A2HOC
                H_[4,2] = H_[2,4]
                H_[3,4] = A2GCBV/EV_[3]
                H_[4,3] = H_[3,4]
                H_[4,4] = -pbm.elpar[iel_][1]*A2HIC+pbm.elpar[iel_][2]*A2HOC
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

