function QCNEW(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : QCNEW
#    *********
# 
#    Source: Quality Control problem 104 from
#    Betty Schultz and Ben Reiser.
# 
#    SIF input: Andrew Conn, August 1992.
#               correction by S. Gratton & Ph. Toint, May 2024
# 
#    classification = "C-OLR2-MY-9-3"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "QCNEW"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["M"] = 9
        v_["1"] = 1
        v_["2"] = 2
        v_["3"] = 3
        v_["4"] = 4
        v_["5"] = 5
        v_["6"] = 6
        v_["7"] = 7
        v_["8"] = 8
        v_["9"] = 9
        v_["10"] = 10
        v_["11"] = 11
        v_["12"] = 12
        v_["13"] = 13
        v_["14"] = 14
        v_["15"] = 15
        v_["16"] = 16
        v_["17"] = 17
        v_["N"] = 10048
        v_["NGG"] = 9900
        v_["NGB"] = 35
        v_["NBG"] = 15
        v_["NBU"] = 18
        v_["NBB"] = 2
        v_["NUU"] = 40
        v_["NUB"] = 21
        v_["NGBUU"] = 9
        v_["NBGUU"] = 2
        v_["NBGBU"] = 2
        v_["NGBUB"] = 3
        v_["NBBUU"] = 1
        v_["NBBBU"] = 0
        v_["NBBUB"] = 0
        v_["S1"] = v_["NGG"]+v_["NGB"]
        v_["S1"] = v_["S1"]+v_["NGBUU"]
        v_["S1"] = v_["S1"]+v_["NGBUB"]
        v_["S2"] = v_["NBG"]+v_["NBB"]
        v_["S2"] = v_["S2"]+v_["NBU"]
        v_["S2"] = v_["S2"]+v_["NBGUU"]
        v_["S2"] = v_["S2"]+v_["NBBUU"]
        v_["S2"] = v_["S2"]+v_["NBGBU"]
        v_["S2"] = v_["S2"]+v_["NBBBU"]
        v_["S2"] = v_["S2"]+v_["NBBUB"]
        v_["S3"] = v_["NGG"]+v_["NBG"]
        v_["S3"] = v_["S3"]+v_["NBGUU"]
        v_["S3"] = v_["S3"]+v_["NBGBU"]
        v_["S4"] = v_["NGB"]+v_["NBB"]
        v_["S4"] = v_["S4"]+v_["NUB"]
        v_["S4"] = v_["S4"]+v_["NGBUU"]
        v_["S4"] = v_["S4"]+v_["NBBUU"]
        v_["S4"] = v_["S4"]+v_["NBBBU"]
        v_["S4"] = v_["S4"]+v_["NGBUB"]
        v_["S4"] = v_["S4"]+v_["NBBUB"]
        v_["U1"] = v_["NGB"]+v_["NGBUU"]
        v_["U1"] = v_["U1"]+v_["NGBUB"]
        v_["L1"] = v_["U1"]+v_["NUU"]
        v_["L1"] = v_["U1"]+v_["NUB"]
        v_["U2"] = v_["NBG"]+v_["NBGUU"]
        v_["U2"] = v_["U2"]+v_["NBGBU"]
        v_["L2"] = v_["U2"]+v_["NUU"]
        v_["L2"] = v_["U2"]+v_["NBU"]
        v_["U3"] = v_["NBB"]+v_["NBBUU"]
        v_["U3"] = v_["U3"]+v_["NBBBU"]
        v_["U3"] = v_["U3"]+v_["NBBUB"]
        v_["L3"] = v_["U3"]+v_["NUU"]
        v_["L3"] = v_["U3"]+v_["NBU"]
        v_["L3"] = v_["U3"]+v_["NUB"]
        v_["L4"] = v_["U2"]+v_["NBBUU"]
        v_["L4"] = v_["L4"]+v_["NBBBU"]
        v_["U4"] = v_["U1"]+v_["NBBUB"]
        v_["U4"] = v_["U4"]+v_["NBU"]
        v_["U4"] = v_["U4"]+v_["NBB"]
        v_["U4"] = v_["U4"]+v_["NUU"]
        v_["U4"] = v_["U4"]+v_["NUB"]
        v_["L5"] = v_["U1"]+v_["NBBUU"]
        v_["L5"] = v_["L5"]+v_["NBBBU"]
        v_["L5"] = v_["L5"]+v_["NBBUB"]
        v_["L5"] = v_["L5"]+v_["NUB"]
        v_["L5"] = v_["L5"]+v_["NBB"]
        v_["U5"] = v_["L5"]+v_["NUU"]
        v_["U5"] = v_["U5"]+v_["NBU"]
        v_["TWO"] = 2.0
        v_["RS1"] = Float64(v_["S1"])
        v_["RS2"] = Float64(v_["S2"])
        v_["RS3"] = Float64(v_["S3"])
        v_["RS4"] = Float64(v_["S4"])
        v_["RNBG"] = Float64(v_["NBG"])
        v_["RNBGBU"] = Float64(v_["NBGBU"])
        v_["RNBGUU"] = Float64(v_["NBGUU"])
        v_["RNGB"] = Float64(v_["NGB"])
        v_["RNGBUB"] = Float64(v_["NGBUB"])
        v_["RNGBUU"] = Float64(v_["NGBUU"])
        v_["RNBB"] = Float64(v_["NBB"])
        v_["RNBBUB"] = Float64(v_["NBBUB"])
        v_["RNBBBU"] = Float64(v_["NBBBU"])
        v_["RNBBUU"] = Float64(v_["NBBUU"])
        v_["RNBGBU"] = Float64(v_["NBGBU"])
        v_["RNUU"] = Float64(v_["NUU"])
        v_["RNUB"] = Float64(v_["NUB"])
        v_["RNBU"] = Float64(v_["NBU"])
        v_["RN"] = Float64(v_["N"])
        v_["RL1"] = Float64(v_["L1"])
        v_["RU1"] = Float64(v_["U1"])
        v_["RL2"] = Float64(v_["L2"])
        v_["RU2"] = Float64(v_["U2"])
        v_["RL3"] = Float64(v_["L3"])
        v_["RU3"] = Float64(v_["U3"])
        v_["RL4"] = Float64(v_["L4"])
        v_["RU4"] = Float64(v_["U4"])
        v_["RL5"] = Float64(v_["L5"])
        v_["RU5"] = Float64(v_["U5"])
        v_["LF1"] = v_["RU4"]/v_["RN"]
        v_["UF1"] = v_["RL4"]/v_["RN"]
        v_["SF1"] = v_["LF1"]+v_["UF1"]
        v_["SF1"] = v_["SF1"]/v_["TWO"]
        v_["LF2"] = v_["RL5"]/v_["RN"]
        v_["UF2"] = v_["RU5"]/v_["RN"]
        v_["SF2"] = v_["LF2"]+v_["UF2"]
        v_["SF2"] = v_["SF2"]/v_["TWO"]
        v_["LGBGB"] = v_["RNGB"]/v_["RL1"]
        v_["UGBGB"] = v_["RNGB"]/v_["RU1"]
        v_["SGBGB"] = v_["LGBGB"]+v_["UGBGB"]
        v_["SGBGB"] = v_["SGBGB"]/v_["TWO"]
        v_["LUBGB"] = v_["RNUB"]/v_["RL1"]
        v_["UUBGB"] = v_["RNUB"]/v_["RU1"]
        v_["SUBGB"] = v_["LUBGB"]+v_["UUBGB"]
        v_["SUBGB"] = v_["SUBGB"]/v_["TWO"]
        v_["LBGBG"] = v_["RNBG"]/v_["RL2"]
        v_["UBGBG"] = v_["RNBG"]/v_["RU2"]
        v_["SBGBG"] = v_["LBGBG"]+v_["UBGBG"]
        v_["SBGBG"] = v_["SBGBG"]/v_["TWO"]
        v_["LBUBG"] = v_["RNBU"]/v_["RL2"]
        v_["UBUBG"] = v_["RNBU"]/v_["RU2"]
        v_["SBUBG"] = v_["LBUBG"]+v_["UBUBG"]
        v_["SBUBG"] = v_["SBUBG"]/v_["TWO"]
        v_["LBBBB"] = v_["RNBB"]/v_["RL3"]
        v_["UBBBB"] = v_["RNBB"]/v_["RU3"]
        v_["SBBBB"] = v_["LBBBB"]+v_["UBBBB"]
        v_["SBBBB"] = v_["SBBBB"]/v_["TWO"]
        v_["LBUBB"] = v_["RNBU"]/v_["RL3"]
        v_["UBUBB"] = v_["RNBU"]/v_["RU3"]
        v_["SBUBB"] = v_["LBUBB"]+v_["UBUBB"]
        v_["SBUBB"] = v_["SBUBB"]/v_["TWO"]
        v_["LUBBB"] = v_["RNUB"]/v_["RL3"]
        v_["UUBBB"] = v_["RNUB"]/v_["RU3"]
        v_["SUBBB"] = v_["LUBBB"]+v_["UUBBB"]
        v_["SUBBB"] = v_["SUBBB"]/v_["TWO"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        iv,ix_,_ = s2mpj_ii("F1",ix_)
        arrset(pb.xnames,iv,"F1")
        iv,ix_,_ = s2mpj_ii("F2",ix_)
        arrset(pb.xnames,iv,"F2")
        iv,ix_,_ = s2mpj_ii("PBGBG",ix_)
        arrset(pb.xnames,iv,"PBGBG")
        iv,ix_,_ = s2mpj_ii("PBUBG",ix_)
        arrset(pb.xnames,iv,"PBUBG")
        iv,ix_,_ = s2mpj_ii("PGBGB",ix_)
        arrset(pb.xnames,iv,"PGBGB")
        iv,ix_,_ = s2mpj_ii("PUBGB",ix_)
        arrset(pb.xnames,iv,"PUBGB")
        iv,ix_,_ = s2mpj_ii("PBBBB",ix_)
        arrset(pb.xnames,iv,"PBBBB")
        iv,ix_,_ = s2mpj_ii("PUBBB",ix_)
        arrset(pb.xnames,iv,"PUBBB")
        iv,ix_,_ = s2mpj_ii("PBUBB",ix_)
        arrset(pb.xnames,iv,"PBUBB")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ"*string(Int64(v_["1"])),ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["F1"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("OBJ"*string(Int64(v_["2"])),ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["F1"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("OBJ"*string(Int64(v_["3"])),ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["F2"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("OBJ"*string(Int64(v_["4"])),ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["F2"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("OBJ"*string(Int64(v_["5"])),ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["PBUBG"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("OBJ"*string(Int64(v_["6"])),ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["PUBGB"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("OBJ"*string(Int64(v_["7"])),ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["F1"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["F2"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("OBJ"*string(Int64(v_["8"])),ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["PGBGB"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("OBJ"*string(Int64(v_["9"])),ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["PBGBG"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("OBJ"*string(Int64(v_["10"])),ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["PBBBB"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("OBJ"*string(Int64(v_["11"])),ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["PGBGB"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["PUBGB"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("OBJ"*string(Int64(v_["12"])),ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["PBGBG"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["PBUBG"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("OBJ"*string(Int64(v_["13"])),ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["PBBBB"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["PUBBB"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("OBJ"*string(Int64(v_["13"])),ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["PBUBB"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("OBJ"*string(Int64(v_["14"])),ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["PBUBG"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("OBJ"*string(Int64(v_["15"])),ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["PBUBB"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("OBJ"*string(Int64(v_["16"])),ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["PUBGB"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("OBJ"*string(Int64(v_["17"])),ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["PUBBB"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("CON0",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CON0")
        iv = ix_["PBGBG"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["PBUBG"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("CON1",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CON1")
        iv = ix_["PGBGB"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["PUBGB"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("CON2",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CON2")
        iv = ix_["PBBBB"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["PUBBB"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["PBUBB"]
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
        pbm.gconst[ig_["OBJ"*string(Int64(v_["1"]))]] = Float64(-1.0)
        pbm.gconst[ig_["OBJ"*string(Int64(v_["3"]))]] = Float64(-1.0)
        pbm.gconst[ig_["OBJ"*string(Int64(v_["7"]))]] = Float64(-1.0)
        pbm.gconst[ig_["OBJ"*string(Int64(v_["11"]))]] = Float64(-1.0)
        pbm.gconst[ig_["OBJ"*string(Int64(v_["12"]))]] = Float64(-1.0)
        pbm.gconst[ig_["OBJ"*string(Int64(v_["13"]))]] = Float64(-1.0)
        pbm.gconst[ig_["CON0"]] = Float64(1.0)
        pbm.gconst[ig_["CON1"]] = Float64(1.0)
        pbm.gconst[ig_["CON2"]] = Float64(1.0)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xupper = fill(1.0,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper[ix_["F1"]] = 8.1608e-3
        pb.xlower[ix_["F1"]] = 7.9618e-3
        pb.xupper[ix_["F2"]] = 1.0649e-2
        pb.xlower[ix_["F2"]] = 1.0450e-2
        pb.xlower[ix_["PBGBG"]] = v_["LBGBG"]
        pb.xupper[ix_["PBGBG"]] = v_["UBGBG"]
        pb.xlower[ix_["PBUBG"]] = v_["LBUBG"]
        pb.xupper[ix_["PBUBG"]] = v_["UBUBG"]
        pb.xlower[ix_["PGBGB"]] = v_["LGBGB"]
        pb.xupper[ix_["PGBGB"]] = v_["UGBGB"]
        pb.xlower[ix_["PUBGB"]] = v_["LUBGB"]
        pb.xupper[ix_["PUBGB"]] = v_["UUBGB"]
        pb.xlower[ix_["PBBBB"]] = v_["LBBBB"]
        pb.xupper[ix_["PBBBB"]] = v_["UBBBB"]
        pb.xlower[ix_["PBUBB"]] = 0.0
        pb.xupper[ix_["PBUBB"]] = 0.0
        pb.xlower[ix_["PUBBB"]] = 0.0
        pb.xupper[ix_["PUBBB"]] = 0.0
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(0.0),pb.n)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "en2PROD", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        it,iet_,_ = s2mpj_ii( "eI2PROD", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        loaset(elftv,it,3,"V3")
        it,iet_,_ = s2mpj_ii( "eI3PROD", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        loaset(elftv,it,3,"V3")
        loaset(elftv,it,4,"V4")
        it,iet_,_ = s2mpj_ii( "en3PRODI", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        loaset(elftv,it,3,"V3")
        loaset(elftv,it,4,"V4")
        loaset(elftv,it,5,"V5")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "E"*string(Int64(v_["1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en2PROD")
        arrset(ielftype,ie,iet_["en2PROD"])
        ename = "E"*string(Int64(v_["1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "F1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.0))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "PBUBG"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.0))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["2"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en2PROD")
        arrset(ielftype,ie,iet_["en2PROD"])
        ename = "E"*string(Int64(v_["2"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "F2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.0))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["2"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "PBUBB"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.0))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["3"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en2PROD")
        arrset(ielftype,ie,iet_["en2PROD"])
        ename = "E"*string(Int64(v_["3"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "F1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.0))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["3"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "PUBGB"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.0))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["4"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en2PROD")
        arrset(ielftype,ie,iet_["en2PROD"])
        ename = "E"*string(Int64(v_["4"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "F1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.0))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["4"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "PUBBB"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.0))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["5"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eI2PROD")
        arrset(ielftype,ie,iet_["eI2PROD"])
        ename = "E"*string(Int64(v_["5"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "F1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.0))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["5"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "PBGBG"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.0))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["5"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "PBUBG"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.0))
        posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["6"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eI3PROD")
        arrset(ielftype,ie,iet_["eI3PROD"])
        ename = "E"*string(Int64(v_["6"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "F1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.0))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["6"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "F2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.0))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["6"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "PBGBG"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.0))
        posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["6"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "PBUBG"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.0))
        posev = findfirst(x->x=="V4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["7"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eI2PROD")
        arrset(ielftype,ie,iet_["eI2PROD"])
        ename = "E"*string(Int64(v_["7"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "F2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.0))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["7"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "PGBGB"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.0))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["7"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "PUBGB"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.0))
        posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["8"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eI3PROD")
        arrset(ielftype,ie,iet_["eI3PROD"])
        ename = "E"*string(Int64(v_["8"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "F1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.0))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["8"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "F2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.0))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["8"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "PGBGB"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.0))
        posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["8"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "PUBGB"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.0))
        posev = findfirst(x->x=="V4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["9"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en3PRODI")
        arrset(ielftype,ie,iet_["en3PRODI"])
        ename = "E"*string(Int64(v_["9"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "F1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.0))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["9"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "F2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.0))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["9"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "PBBBB"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.0))
        posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["9"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "PUBBB"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.0))
        posev = findfirst(x->x=="V4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["9"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "PBUBB"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.0))
        posev = findfirst(x->x=="V5",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gLOG",igt_)
        it,igt_,_ = s2mpj_ii("gLOG",igt_)
        grftp = Vector{Vector{String}}()
        loaset(grftp,it,1,"P")
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["OBJ"*string(Int64(v_["1"]))]
        arrset(pbm.grftype,ig,"gLOG")
        posgp = findfirst(x->x=="P",grftp[igt_[pbm.grftype[ig]]])
        loaset(pbm.grpar,ig,posgp,Float64(v_["RS1"]))
        ig = ig_["OBJ"*string(Int64(v_["2"]))]
        arrset(pbm.grftype,ig,"gLOG")
        posgp = findfirst(x->x=="P",grftp[igt_[pbm.grftype[ig]]])
        loaset(pbm.grpar,ig,posgp,Float64(v_["RS2"]))
        ig = ig_["OBJ"*string(Int64(v_["3"]))]
        arrset(pbm.grftype,ig,"gLOG")
        posgp = findfirst(x->x=="P",grftp[igt_[pbm.grftype[ig]]])
        loaset(pbm.grpar,ig,posgp,Float64(v_["RS3"]))
        ig = ig_["OBJ"*string(Int64(v_["4"]))]
        arrset(pbm.grftype,ig,"gLOG")
        posgp = findfirst(x->x=="P",grftp[igt_[pbm.grftype[ig]]])
        loaset(pbm.grpar,ig,posgp,Float64(v_["RS4"]))
        ig = ig_["OBJ"*string(Int64(v_["5"]))]
        arrset(pbm.grftype,ig,"gLOG")
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["1"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["2"]))])
        loaset(pbm.grelw,ig,posel, 1.)
        posgp = findfirst(x->x=="P",grftp[igt_[pbm.grftype[ig]]])
        loaset(pbm.grpar,ig,posgp,Float64(v_["RNBU"]))
        ig = ig_["OBJ"*string(Int64(v_["6"]))]
        arrset(pbm.grftype,ig,"gLOG")
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["3"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["4"]))])
        loaset(pbm.grelw,ig,posel, 1.)
        posgp = findfirst(x->x=="P",grftp[igt_[pbm.grftype[ig]]])
        loaset(pbm.grpar,ig,posgp,Float64(v_["RNUB"]))
        ig = ig_["OBJ"*string(Int64(v_["7"]))]
        arrset(pbm.grftype,ig,"gLOG")
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["5"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["6"]))])
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["7"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["8"]))])
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["9"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posgp = findfirst(x->x=="P",grftp[igt_[pbm.grftype[ig]]])
        loaset(pbm.grpar,ig,posgp,Float64(v_["RNUU"]))
        ig = ig_["OBJ"*string(Int64(v_["8"]))]
        arrset(pbm.grftype,ig,"gLOG")
        posgp = findfirst(x->x=="P",grftp[igt_[pbm.grftype[ig]]])
        loaset(pbm.grpar,ig,posgp,Float64(v_["RNGB"]))
        ig = ig_["OBJ"*string(Int64(v_["9"]))]
        arrset(pbm.grftype,ig,"gLOG")
        posgp = findfirst(x->x=="P",grftp[igt_[pbm.grftype[ig]]])
        loaset(pbm.grpar,ig,posgp,Float64(v_["RNBG"]))
        ig = ig_["OBJ"*string(Int64(v_["10"]))]
        arrset(pbm.grftype,ig,"gLOG")
        posgp = findfirst(x->x=="P",grftp[igt_[pbm.grftype[ig]]])
        loaset(pbm.grpar,ig,posgp,Float64(v_["RNBB"]))
        ig = ig_["OBJ"*string(Int64(v_["11"]))]
        arrset(pbm.grftype,ig,"gLOG")
        posgp = findfirst(x->x=="P",grftp[igt_[pbm.grftype[ig]]])
        loaset(pbm.grpar,ig,posgp,Float64(v_["RNGBUU"]))
        ig = ig_["OBJ"*string(Int64(v_["12"]))]
        arrset(pbm.grftype,ig,"gLOG")
        posgp = findfirst(x->x=="P",grftp[igt_[pbm.grftype[ig]]])
        loaset(pbm.grpar,ig,posgp,Float64(v_["RNBGUU"]))
        ig = ig_["OBJ"*string(Int64(v_["13"]))]
        arrset(pbm.grftype,ig,"gLOG")
        posgp = findfirst(x->x=="P",grftp[igt_[pbm.grftype[ig]]])
        loaset(pbm.grpar,ig,posgp,Float64(v_["RNBBUU"]))
        ig = ig_["OBJ"*string(Int64(v_["14"]))]
        arrset(pbm.grftype,ig,"gLOG")
        posgp = findfirst(x->x=="P",grftp[igt_[pbm.grftype[ig]]])
        loaset(pbm.grpar,ig,posgp,Float64(v_["RNBGBU"]))
        ig = ig_["OBJ"*string(Int64(v_["15"]))]
        arrset(pbm.grftype,ig,"gLOG")
        posgp = findfirst(x->x=="P",grftp[igt_[pbm.grftype[ig]]])
        loaset(pbm.grpar,ig,posgp,Float64(v_["RNBBBU"]))
        ig = ig_["OBJ"*string(Int64(v_["16"]))]
        arrset(pbm.grftype,ig,"gLOG")
        posgp = findfirst(x->x=="P",grftp[igt_[pbm.grftype[ig]]])
        loaset(pbm.grpar,ig,posgp,Float64(v_["RNGBUB"]))
        ig = ig_["OBJ"*string(Int64(v_["17"]))]
        arrset(pbm.grftype,ig,"gLOG")
        posgp = findfirst(x->x=="P",grftp[igt_[pbm.grftype[ig]]])
        loaset(pbm.grpar,ig,posgp,Float64(v_["RNBBUU"]))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               1138.416240
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-OLR2-MY-9-3"
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

    elseif action == "en2PROD"

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
                H_[1,1] = 0.0
                H_[2,2] = 0.0
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

    elseif action == "eI2PROD"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,2,3)
        IV_ =  zeros(Float64,2)
        U_[1,1] = U_[1,1]+1
        U_[2,2] = U_[2,2]-1
        U_[2,3] = U_[2,3]-1
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        f_   = IV_[1]*(1.0+IV_[2])
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 1.0+IV_[2]
            g_[2] = IV_[1]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = 0.0
                H_[2,2] = 0.0
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

    elseif action == "eI3PROD"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,3,4)
        IV_ =  zeros(Float64,3)
        U_[1,1] = U_[1,1]+1
        U_[2,2] = U_[2,2]+1
        U_[3,3] = U_[3,3]-1
        U_[3,4] = U_[3,4]-1
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        IV_[3] = dot(U_[3,:],EV_)
        f_   = IV_[1]*IV_[2]*(1.0+IV_[3])
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = IV_[2]*(1.0+IV_[3])
            g_[2] = IV_[1]*(1.0+IV_[3])
            g_[3] = IV_[1]*IV_[2]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,1] = 0.0
                H_[1,2] = 1.0+IV_[3]
                H_[2,1] = H_[1,2]
                H_[1,3] = IV_[2]
                H_[3,1] = H_[1,3]
                H_[2,2] = 0.0
                H_[2,3] = IV_[1]
                H_[3,2] = H_[2,3]
                H_[3,3] = 0.0
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

    elseif action == "en3PRODI"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,3,5)
        IV_ =  zeros(Float64,3)
        U_[1,1] = U_[1,1]+1
        U_[2,2] = U_[2,2]+1
        U_[3,3] = U_[3,3]-1
        U_[3,4] = U_[3,4]-1
        U_[3,5] = U_[3,5]-1
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        IV_[3] = dot(U_[3,:],EV_)
        f_   = IV_[1]*IV_[2]*(1.0+IV_[3])
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = IV_[2]*(1.0+IV_[3])
            g_[2] = IV_[1]*(1.0+IV_[3])
            g_[3] = IV_[1]*IV_[2]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,1] = 0.0
                H_[1,2] = 1.0+IV_[3]
                H_[2,1] = H_[1,2]
                H_[1,3] = IV_[2]
                H_[3,1] = H_[1,3]
                H_[2,2] = 0.0
                H_[2,3] = IV_[1]
                H_[3,2] = H_[2,3]
                H_[3,3] = 0.0
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

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    elseif action == "gLOG"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        ARG0 = GVAR_==0.0
        if ARG0
            FF = 0.0
        end
        if !ARG0
            FF = pbm.grpar[igr_][1]*log(GVAR_)
        end
        if ARG0
            GG = 0.0
        end
        if !ARG0
            GG = pbm.grpar[igr_][1]/GVAR_
        end
        if ARG0
            HH = 0.0
        end
        if !ARG0
            HH = -pbm.grpar[igr_][1]/GVAR_^2
        end
        f_= FF
        if nargout>1
            g_ = GG
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = HH
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

