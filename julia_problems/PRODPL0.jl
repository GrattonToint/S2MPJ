function PRODPL0(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : PRODPL0
#    *********
# 
#    A production planning problem in the computer industry.
# 
#    Source:
#    L. Escudero, private communication, 1991.
# 
#    SIF input: A.R. Conn, March 1991.
# 
#    classification = "C-LQR2-RY-60-29"
# 
#    Constants
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "PRODPL0"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["T"] = 5
        v_["1"] = 1
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("COST",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("K01",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"K01")
        ig,ig_,_ = s2mpj_ii("K02",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"K02")
        ig,ig_,_ = s2mpj_ii("K03",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"K03")
        ig,ig_,_ = s2mpj_ii("K04",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"K04")
        ig,ig_,_ = s2mpj_ii("K05",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"K05")
        ig,ig_,_ = s2mpj_ii("D00101",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"D00101")
        ig,ig_,_ = s2mpj_ii("D00201",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"D00201")
        ig,ig_,_ = s2mpj_ii("D00301",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"D00301")
        ig,ig_,_ = s2mpj_ii("D00401",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"D00401")
        ig,ig_,_ = s2mpj_ii("D00102",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"D00102")
        ig,ig_,_ = s2mpj_ii("D00202",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"D00202")
        ig,ig_,_ = s2mpj_ii("D00302",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"D00302")
        ig,ig_,_ = s2mpj_ii("D00402",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"D00402")
        ig,ig_,_ = s2mpj_ii("D00103",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"D00103")
        ig,ig_,_ = s2mpj_ii("D00203",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"D00203")
        ig,ig_,_ = s2mpj_ii("D00303",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"D00303")
        ig,ig_,_ = s2mpj_ii("D00403",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"D00403")
        ig,ig_,_ = s2mpj_ii("D00104",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"D00104")
        ig,ig_,_ = s2mpj_ii("D00204",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"D00204")
        ig,ig_,_ = s2mpj_ii("D00304",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"D00304")
        ig,ig_,_ = s2mpj_ii("D00404",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"D00404")
        ig,ig_,_ = s2mpj_ii("D00105",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"D00105")
        ig,ig_,_ = s2mpj_ii("D00205",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"D00205")
        ig,ig_,_ = s2mpj_ii("D00305",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"D00305")
        ig,ig_,_ = s2mpj_ii("D00405",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"D00405")
        v_["TM1"] = -1+v_["T"]
        for I = Int64(v_["1"]):Int64(v_["TM1"])
            ig,ig_,_ = s2mpj_ii("SMOOTH"*string(I),ig_)
            arrset(gtype,ig,"<=")
            arrset(pb.cnames,ig,"SMOOTH"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        ngrp   = length(ig_)
        iv,ix_,_ = s2mpj_ii("X00101",ix_)
        arrset(pb.xnames,iv,"X00101")
        ig = ig_["K01"]
        pbm.A[ig,iv] += Float64(1)
        ig = ig_["D00101"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X00201",ix_)
        arrset(pb.xnames,iv,"X00201")
        ig = ig_["K01"]
        pbm.A[ig,iv] += Float64(1)
        ig = ig_["D00201"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X00301",ix_)
        arrset(pb.xnames,iv,"X00301")
        ig = ig_["K01"]
        pbm.A[ig,iv] += Float64(1)
        ig = ig_["D00301"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X00401",ix_)
        arrset(pb.xnames,iv,"X00401")
        ig = ig_["K01"]
        pbm.A[ig,iv] += Float64(1)
        ig = ig_["D00401"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X00102",ix_)
        arrset(pb.xnames,iv,"X00102")
        ig = ig_["K02"]
        pbm.A[ig,iv] += Float64(1)
        ig = ig_["D00102"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X00202",ix_)
        arrset(pb.xnames,iv,"X00202")
        ig = ig_["K02"]
        pbm.A[ig,iv] += Float64(1)
        ig = ig_["D00202"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X00302",ix_)
        arrset(pb.xnames,iv,"X00302")
        ig = ig_["K02"]
        pbm.A[ig,iv] += Float64(1)
        ig = ig_["D00302"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X00402",ix_)
        arrset(pb.xnames,iv,"X00402")
        ig = ig_["K02"]
        pbm.A[ig,iv] += Float64(1)
        ig = ig_["D00402"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X00103",ix_)
        arrset(pb.xnames,iv,"X00103")
        ig = ig_["K03"]
        pbm.A[ig,iv] += Float64(1)
        ig = ig_["D00103"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X00203",ix_)
        arrset(pb.xnames,iv,"X00203")
        ig = ig_["K03"]
        pbm.A[ig,iv] += Float64(1)
        ig = ig_["D00203"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X00303",ix_)
        arrset(pb.xnames,iv,"X00303")
        ig = ig_["K03"]
        pbm.A[ig,iv] += Float64(1)
        ig = ig_["D00303"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X00403",ix_)
        arrset(pb.xnames,iv,"X00403")
        ig = ig_["K03"]
        pbm.A[ig,iv] += Float64(1)
        ig = ig_["D00403"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X00104",ix_)
        arrset(pb.xnames,iv,"X00104")
        ig = ig_["K04"]
        pbm.A[ig,iv] += Float64(1)
        ig = ig_["D00104"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X00204",ix_)
        arrset(pb.xnames,iv,"X00204")
        ig = ig_["K04"]
        pbm.A[ig,iv] += Float64(1)
        ig = ig_["D00204"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X00304",ix_)
        arrset(pb.xnames,iv,"X00304")
        ig = ig_["K04"]
        pbm.A[ig,iv] += Float64(1)
        ig = ig_["D00304"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X00404",ix_)
        arrset(pb.xnames,iv,"X00404")
        ig = ig_["K04"]
        pbm.A[ig,iv] += Float64(1)
        ig = ig_["D00404"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X00105",ix_)
        arrset(pb.xnames,iv,"X00105")
        ig = ig_["K05"]
        pbm.A[ig,iv] += Float64(1)
        ig = ig_["D00105"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X00205",ix_)
        arrset(pb.xnames,iv,"X00205")
        ig = ig_["K05"]
        pbm.A[ig,iv] += Float64(1)
        ig = ig_["D00205"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X00305",ix_)
        arrset(pb.xnames,iv,"X00305")
        ig = ig_["K05"]
        pbm.A[ig,iv] += Float64(1)
        ig = ig_["D00305"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("X00405",ix_)
        arrset(pb.xnames,iv,"X00405")
        ig = ig_["K05"]
        pbm.A[ig,iv] += Float64(1)
        ig = ig_["D00405"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00101",ix_)
        arrset(pb.xnames,iv,"I00101")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(1.000000)
        ig = ig_["D00101"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00101",ix_)
        arrset(pb.xnames,iv,"I00101")
        ig = ig_["D00102"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00101",ix_)
        arrset(pb.xnames,iv,"Y00101")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(2)
        ig = ig_["D00101"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00201",ix_)
        arrset(pb.xnames,iv,"I00201")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(2.000000)
        ig = ig_["D00201"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00201",ix_)
        arrset(pb.xnames,iv,"I00201")
        ig = ig_["D00202"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00201",ix_)
        arrset(pb.xnames,iv,"Y00201")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(3)
        ig = ig_["D00201"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00301",ix_)
        arrset(pb.xnames,iv,"I00301")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(3.000000)
        ig = ig_["D00301"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00301",ix_)
        arrset(pb.xnames,iv,"I00301")
        ig = ig_["D00302"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00301",ix_)
        arrset(pb.xnames,iv,"Y00301")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(2)
        ig = ig_["D00301"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00401",ix_)
        arrset(pb.xnames,iv,"I00401")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(4.000000)
        ig = ig_["D00401"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00401",ix_)
        arrset(pb.xnames,iv,"I00401")
        ig = ig_["D00402"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00401",ix_)
        arrset(pb.xnames,iv,"Y00401")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(5)
        ig = ig_["D00401"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00102",ix_)
        arrset(pb.xnames,iv,"I00102")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(1.000000)
        ig = ig_["D00102"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00102",ix_)
        arrset(pb.xnames,iv,"I00102")
        ig = ig_["D00103"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00102",ix_)
        arrset(pb.xnames,iv,"Y00102")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(2)
        ig = ig_["D00102"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00202",ix_)
        arrset(pb.xnames,iv,"I00202")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(2.000000)
        ig = ig_["D00202"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00202",ix_)
        arrset(pb.xnames,iv,"I00202")
        ig = ig_["D00203"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00202",ix_)
        arrset(pb.xnames,iv,"Y00202")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(3)
        ig = ig_["D00202"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00302",ix_)
        arrset(pb.xnames,iv,"I00302")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(3.000000)
        ig = ig_["D00302"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00302",ix_)
        arrset(pb.xnames,iv,"I00302")
        ig = ig_["D00303"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00302",ix_)
        arrset(pb.xnames,iv,"Y00302")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(2)
        ig = ig_["D00302"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00402",ix_)
        arrset(pb.xnames,iv,"I00402")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(4.000000)
        ig = ig_["D00402"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00402",ix_)
        arrset(pb.xnames,iv,"I00402")
        ig = ig_["D00403"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00402",ix_)
        arrset(pb.xnames,iv,"Y00402")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(5)
        ig = ig_["D00402"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00103",ix_)
        arrset(pb.xnames,iv,"I00103")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(1.000000)
        ig = ig_["D00103"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00103",ix_)
        arrset(pb.xnames,iv,"I00103")
        ig = ig_["D00104"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00103",ix_)
        arrset(pb.xnames,iv,"Y00103")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(2)
        ig = ig_["D00103"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00203",ix_)
        arrset(pb.xnames,iv,"I00203")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(2.000000)
        ig = ig_["D00203"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00203",ix_)
        arrset(pb.xnames,iv,"I00203")
        ig = ig_["D00204"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00203",ix_)
        arrset(pb.xnames,iv,"Y00203")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(3)
        ig = ig_["D00203"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00303",ix_)
        arrset(pb.xnames,iv,"I00303")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(3.000000)
        ig = ig_["D00303"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00303",ix_)
        arrset(pb.xnames,iv,"I00303")
        ig = ig_["D00304"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00303",ix_)
        arrset(pb.xnames,iv,"Y00303")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(2)
        ig = ig_["D00303"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00403",ix_)
        arrset(pb.xnames,iv,"I00403")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(4.000000)
        ig = ig_["D00403"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00403",ix_)
        arrset(pb.xnames,iv,"I00403")
        ig = ig_["D00404"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00403",ix_)
        arrset(pb.xnames,iv,"Y00403")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(5)
        ig = ig_["D00403"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00104",ix_)
        arrset(pb.xnames,iv,"I00104")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(1.000000)
        ig = ig_["D00104"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00104",ix_)
        arrset(pb.xnames,iv,"I00104")
        ig = ig_["D00105"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00104",ix_)
        arrset(pb.xnames,iv,"Y00104")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(2)
        ig = ig_["D00104"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00204",ix_)
        arrset(pb.xnames,iv,"I00204")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(2.000000)
        ig = ig_["D00204"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00204",ix_)
        arrset(pb.xnames,iv,"I00204")
        ig = ig_["D00205"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00204",ix_)
        arrset(pb.xnames,iv,"Y00204")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(3)
        ig = ig_["D00204"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00304",ix_)
        arrset(pb.xnames,iv,"I00304")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(3.000000)
        ig = ig_["D00304"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00304",ix_)
        arrset(pb.xnames,iv,"I00304")
        ig = ig_["D00305"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00304",ix_)
        arrset(pb.xnames,iv,"Y00304")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(2)
        ig = ig_["D00304"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00404",ix_)
        arrset(pb.xnames,iv,"I00404")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(4.000000)
        ig = ig_["D00404"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("I00404",ix_)
        arrset(pb.xnames,iv,"I00404")
        ig = ig_["D00405"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("Y00404",ix_)
        arrset(pb.xnames,iv,"Y00404")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(5)
        ig = ig_["D00404"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00105",ix_)
        arrset(pb.xnames,iv,"I00105")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(1.000000)
        ig = ig_["D00105"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("Y00105",ix_)
        arrset(pb.xnames,iv,"Y00105")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(2)
        ig = ig_["D00105"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00205",ix_)
        arrset(pb.xnames,iv,"I00205")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(2.000000)
        ig = ig_["D00205"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("Y00205",ix_)
        arrset(pb.xnames,iv,"Y00205")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(3)
        ig = ig_["D00205"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00305",ix_)
        arrset(pb.xnames,iv,"I00305")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(3.000000)
        ig = ig_["D00305"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("Y00305",ix_)
        arrset(pb.xnames,iv,"Y00305")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(2)
        ig = ig_["D00305"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("I00405",ix_)
        arrset(pb.xnames,iv,"I00405")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(4.000000)
        ig = ig_["D00405"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("Y00405",ix_)
        arrset(pb.xnames,iv,"Y00405")
        ig = ig_["COST"]
        pbm.A[ig,iv] += Float64(5)
        ig = ig_["D00405"]
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
        pbm.gconst[ig_["K01"]] = Float64(3)
        pbm.gconst[ig_["K02"]] = Float64(6)
        pbm.gconst[ig_["K03"]] = Float64(10)
        pbm.gconst[ig_["K04"]] = Float64(2000)
        pbm.gconst[ig_["K05"]] = Float64(18)
        pbm.gconst[ig_["D00101"]] = Float64(1.000)
        pbm.gconst[ig_["D00201"]] = Float64(1.000)
        pbm.gconst[ig_["D00301"]] = Float64(1.000)
        pbm.gconst[ig_["D00401"]] = Float64(1.000)
        pbm.gconst[ig_["D00102"]] = Float64(2.667)
        pbm.gconst[ig_["D00202"]] = Float64(1.667)
        pbm.gconst[ig_["D00302"]] = Float64(2.667)
        pbm.gconst[ig_["D00402"]] = Float64(3.333)
        pbm.gconst[ig_["D00103"]] = Float64(2.667)
        pbm.gconst[ig_["D00203"]] = Float64(2.000)
        pbm.gconst[ig_["D00303"]] = Float64(3.000)
        pbm.gconst[ig_["D00403"]] = Float64(3.000)
        pbm.gconst[ig_["D00104"]] = Float64(2.667)
        pbm.gconst[ig_["D00204"]] = Float64(2.667)
        pbm.gconst[ig_["D00304"]] = Float64(2.667)
        pbm.gconst[ig_["D00404"]] = Float64(2.667)
        pbm.gconst[ig_["D00105"]] = Float64(2.667)
        pbm.gconst[ig_["D00205"]] = Float64(2.333)
        pbm.gconst[ig_["D00305"]] = Float64(2.333)
        pbm.gconst[ig_["D00405"]] = Float64(2.333)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSQMRSQ", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        loaset(elftv,it,3,"V3")
        loaset(elftv,it,4,"V4")
        loaset(elftv,it,5,"V5")
        loaset(elftv,it,6,"V6")
        loaset(elftv,it,7,"V7")
        loaset(elftv,it,8,"V8")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["TM1"])
            v_["IP1"] = 1+I
            ename = "NLE"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQMRSQ")
            arrset(ielftype,ie,iet_["eSQMRSQ"])
            vname = "X0010"*string(Int64(v_["IP1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X0020"*string(Int64(v_["IP1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X0030"*string(Int64(v_["IP1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X0040"*string(Int64(v_["IP1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V4",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X0010"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V5",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X0020"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V6",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X0030"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V7",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X0040"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V8",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["TM1"])
            ig = ig_["SMOOTH"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["NLE"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               58.7898356794
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
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
        pb.pbclass = "C-LQR2-RY-60-29"
        pb.x0          = zeros(Float64,pb.n)
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

    elseif action == "e_globs"

        pbm = args[1]
        arrset(pbm.efpar,1,2.0)
        arrset(pbm.efpar,2,0.1)
        arrset(pbm.efpar,3,pbm.efpar[2]*pbm.efpar[2])
        return pbm

    elseif action == "eSQMRSQ"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,2,8)
        IV_ =  zeros(Float64,2)
        U_[1,1] = U_[1,1]+1
        U_[1,2] = U_[1,2]+1
        U_[1,3] = U_[1,3]+1
        U_[1,4] = U_[1,4]+1
        U_[2,5] = U_[2,5]+1
        U_[2,6] = U_[2,6]+1
        U_[2,7] = U_[2,7]+1
        U_[2,8] = U_[2,8]+1
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        U1MU2 = IV_[1]-IV_[2]
        f_   = U1MU2^2-pbm.efpar[3]*IV_[2]^2
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = pbm.efpar[1]*U1MU2
            g_[2] = -pbm.efpar[1]*U1MU2-pbm.efpar[1]*pbm.efpar[3]*IV_[2]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = pbm.efpar[1]
                H_[1,2] = -pbm.efpar[1]
                H_[2,1] = H_[1,2]
                H_[2,2] = pbm.efpar[1]*(1.0-pbm.efpar[3])
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

    #%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    elseif action in  ["fx","fgx","fgHx","cx","cJx","cJHx","cIx","cIJx","cIJHx","cIJxv","fHxv",
                       "cJxv","cJtxv","cIJtxv","Lxy","Lgxy","LgHxy","LIxy","LIgxy","LIgHxy",
                       "LHxyv","LIHxyv"]

        pbm = args[1]
        if pbm.name == name
            pbm.has_globs = [3,0]
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

