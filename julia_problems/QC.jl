function QC(action::String,args::Union{Any}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : QC
#    *********
# 
#    Source: Quality Control problem 104 from
#    Betty Schultz and Ben Reiser.
# 
#    SIF input: Andrew Conn, August 1992.
#               correction by S. Gratton & Ph. Toint, May 2024
# 
#    classification = "C-COLR2-MY-9-4"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 21 VI 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "QC"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling QC.")
    end

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
        v_["MM"] = v_["N"]-v_["NGG"]
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
        v_["L1"] = v_["L1"]+v_["NUB"]
        v_["U2"] = v_["NBG"]+v_["NBGUU"]
        v_["U2"] = v_["U2"]+v_["NBGBU"]
        v_["L2"] = v_["U2"]+v_["NUU"]
        v_["L2"] = v_["L2"]+v_["NBU"]
        v_["U3"] = v_["NBB"]+v_["NBBUU"]
        v_["U3"] = v_["U3"]+v_["NBBBU"]
        v_["U3"] = v_["U3"]+v_["NBBUB"]
        v_["L3"] = v_["U3"]+v_["NUU"]
        v_["L3"] = v_["L3"]+v_["NBU"]
        v_["L3"] = v_["L3"]+v_["NUB"]
        v_["U4"] = v_["U2"]+v_["NBU"]
        v_["L4"] = v_["U4"]+v_["NUU"]
        v_["U5"] = v_["U1"]+v_["NUB"]
        v_["L5"] = v_["U5"]+v_["NUU"]
        v_["U6"] = v_["U3"]+v_["NUB"]
        v_["L6"] = v_["U6"]+v_["NBU"]
        v_["L6"] = v_["L6"]+v_["NUU"]
        v_["U7"] = v_["U3"]+v_["NBU"]
        v_["L7"] = v_["U7"]+v_["NUB"]
        v_["L7"] = v_["L7"]+v_["NUU"]
        v_["L8"] = v_["U2"]+v_["NBBUU"]
        v_["L8"] = v_["L8"]+v_["NBBBU"]
        v_["L8"] = v_["L8"]+v_["NBBUB"]
        v_["L8"] = v_["L8"]+v_["NBU"]
        v_["L8"] = v_["L8"]+v_["NBB"]
        v_["U8"] = v_["L8"]+v_["NUU"]
        v_["U8"] = v_["U8"]+v_["NUB"]
        v_["L9"] = v_["U1"]+v_["NBBUU"]
        v_["L9"] = v_["L9"]+v_["NBBBU"]
        v_["L9"] = v_["L9"]+v_["NBBUB"]
        v_["L9"] = v_["L9"]+v_["NUB"]
        v_["L9"] = v_["L9"]+v_["NBB"]
        v_["U9"] = v_["L9"]+v_["NUU"]
        v_["U9"] = v_["U9"]+v_["NBU"]
        v_["ZERO"] = 0.0
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
        v_["RL6"] = Float64(v_["L6"])
        v_["RU6"] = Float64(v_["U6"])
        v_["RL7"] = Float64(v_["L7"])
        v_["RU7"] = Float64(v_["U7"])
        v_["RL8"] = Float64(v_["L8"])
        v_["RU8"] = Float64(v_["U8"])
        v_["RL9"] = Float64(v_["L9"])
        v_["RU9"] = Float64(v_["U9"])
        v_["LF1"] = v_["RL8"]/v_["RN"]
        v_["UF1"] = v_["RU8"]/v_["RN"]
        v_["SF1"] = v_["LF1"]+v_["UF1"]
        v_["SF1"] = v_["SF1"]/v_["TWO"]
        v_["LF2"] = v_["RL9"]/v_["RN"]
        v_["UF2"] = v_["RU9"]/v_["RN"]
        v_["SF2"] = v_["LF2"]+v_["UF2"]
        v_["SF2"] = v_["SF2"]/v_["TWO"]
        v_["LGBGB"] = v_["RNGB"]/v_["RL1"]
        v_["UGBGB"] = v_["RNGB"]/v_["RU1"]
        v_["SGBGB"] = v_["LGBGB"]+v_["UGBGB"]
        v_["SGBGB"] = v_["SGBGB"]/v_["TWO"]
        v_["LUBGB"] = v_["RNGBUB"]/v_["RL5"]
        v_["UUBGB"] = v_["RNGBUB"]+v_["RNUB"]
        v_["UUBGB"] = v_["UUBGB"]/v_["RU5"]
        v_["SUBGB"] = v_["LUBGB"]+v_["UUBGB"]
        v_["SUBGB"] = v_["SUBGB"]/v_["TWO"]
        v_["LBGBG"] = v_["RNBG"]/v_["RL2"]
        v_["UBGBG"] = v_["RNBG"]/v_["RU2"]
        v_["SBGBG"] = v_["LBGBG"]+v_["UBGBG"]
        v_["SBGBG"] = v_["SBGBG"]/v_["TWO"]
        v_["LBUBG"] = v_["RNBGBU"]/v_["RL4"]
        v_["UBUBG"] = v_["RNBGBU"]+v_["RNBU"]
        v_["UBUBG"] = v_["UBUBG"]/v_["RU4"]
        v_["SBUBG"] = v_["LBUBG"]+v_["UBUBG"]
        v_["SBUBG"] = v_["SBUBG"]/v_["TWO"]
        v_["LBBBB"] = v_["RNBB"]/v_["RL3"]
        v_["UBBBB"] = v_["RNBB"]/v_["RU3"]
        v_["SBBBB"] = v_["LBBBB"]+v_["UBBBB"]
        v_["SBBBB"] = v_["SBBBB"]/v_["TWO"]
        v_["LBUBB"] = v_["RNBBBU"]/v_["RL7"]
        v_["UBUBB"] = v_["RNBBBU"]+v_["RNBU"]
        v_["UBUBB"] = v_["UBUBB"]/v_["RU7"]
        v_["SBUBB"] = v_["LBUBB"]+v_["UBUBB"]
        v_["SBUBB"] = v_["SBUBB"]/v_["TWO"]
        v_["LUBBB"] = v_["RNBBUB"]/v_["RL6"]
        v_["UUBBB"] = v_["RNBBUB"]+v_["RNUB"]
        v_["UUBBB"] = v_["UUBBB"]/v_["RU6"]
        v_["SUBBB"] = v_["LUBBB"]+v_["UUBBB"]
        v_["SUBBB"] = v_["SUBBB"]/v_["TWO"]
        v_["RMM"] = Float64(v_["MM"])
        v_["RMM"] = v_["RMM"]/v_["RN"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        irA   = Int64[]
        icA   = Int64[]
        valA  = Float64[]
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
        gtype = String[]
        ig,ig_,_ = s2mpj_ii("OBJ"*string(Int64(v_["1"])),ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["F1"])
        push!(valA,Float64(-1.0))
        ig,ig_,_ = s2mpj_ii("OBJ"*string(Int64(v_["2"])),ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["F1"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("OBJ"*string(Int64(v_["3"])),ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["F2"])
        push!(valA,Float64(-1.0))
        ig,ig_,_ = s2mpj_ii("OBJ"*string(Int64(v_["4"])),ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["F2"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("OBJ"*string(Int64(v_["5"])),ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["PBUBG"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("OBJ"*string(Int64(v_["6"])),ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["PUBGB"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("OBJ"*string(Int64(v_["7"])),ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["F1"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["F2"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("OBJ"*string(Int64(v_["8"])),ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["PGBGB"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("OBJ"*string(Int64(v_["9"])),ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["PBGBG"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("OBJ"*string(Int64(v_["10"])),ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["PBBBB"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("OBJ"*string(Int64(v_["11"])),ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["PGBGB"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["PUBGB"])
        push!(valA,Float64(-1.0))
        ig,ig_,_ = s2mpj_ii("OBJ"*string(Int64(v_["12"])),ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["PBGBG"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["PBUBG"])
        push!(valA,Float64(-1.0))
        ig,ig_,_ = s2mpj_ii("OBJ"*string(Int64(v_["13"])),ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["PBBBB"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["PUBBB"])
        push!(valA,Float64(-1.0))
        ig,ig_,_ = s2mpj_ii("OBJ"*string(Int64(v_["13"])),ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["PBUBB"])
        push!(valA,Float64(-1.0))
        ig,ig_,_ = s2mpj_ii("OBJ"*string(Int64(v_["14"])),ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["PBUBG"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("OBJ"*string(Int64(v_["15"])),ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["PBUBB"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("OBJ"*string(Int64(v_["16"])),ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["PUBGB"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("OBJ"*string(Int64(v_["17"])),ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["PUBBB"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("CON0",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CON0")
        push!(irA,ig)
        push!(icA,ix_["PBGBG"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["PBUBG"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("CON1",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CON1")
        push!(irA,ig)
        push!(icA,ix_["PGBGB"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["PUBGB"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("CON2",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CON2")
        push!(irA,ig)
        push!(icA,ix_["PBBBB"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["PUBBB"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["PBUBB"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("CON3",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"CON3")
        push!(irA,ig)
        push!(icA,ix_["F1"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["F2"])
        push!(valA,Float64(1.0))
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
        pbm.gconst[ig_["CON3"]] = Float64(v_["RMM"])
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xupper = fill(1.0,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        pb.xlower[ix_["F1"]] = v_["LF1"]
        pb.xupper[ix_["F1"]] = v_["UF1"]
        pb.xlower[ix_["F2"]] = v_["LF2"]
        pb.xupper[ix_["F2"]] = v_["UF2"]
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
        pb.x0 = fill(Float64(1.5),pb.n)
        pb.x0[ix_["F1"]] = Float64(v_["SF1"])
        pb.x0[ix_["F2"]] = Float64(v_["SF2"])
        pb.x0[ix_["PBGBG"]] = Float64(v_["SBGBG"])
        pb.x0[ix_["PBUBG"]] = Float64(v_["SBUBG"])
        pb.x0[ix_["PGBGB"]] = Float64(v_["SGBGB"])
        pb.x0[ix_["PUBGB"]] = Float64(v_["SUBGB"])
        pb.x0[ix_["PBBBB"]] = Float64(v_["SBBBB"])
        pb.x0[ix_["PBUBB"]] = Float64(v_["ZERO"])
        pb.x0[ix_["PUBBB"]] = Float64(v_["ZERO"])
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
        vname = "F2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(1.5))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "PBUBG"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(1.5))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["2"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en2PROD")
        arrset(ielftype,ie,iet_["en2PROD"])
        ename = "E"*string(Int64(v_["2"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "F2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(1.5))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["2"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "PBUBB"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(1.5))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["3"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en2PROD")
        arrset(ielftype,ie,iet_["en2PROD"])
        ename = "E"*string(Int64(v_["3"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "F1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(1.5))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["3"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "PUBGB"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(1.5))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["4"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en2PROD")
        arrset(ielftype,ie,iet_["en2PROD"])
        ename = "E"*string(Int64(v_["4"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "F1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(1.5))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["4"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "PUBBB"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(1.5))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["5"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eI2PROD")
        arrset(ielftype,ie,iet_["eI2PROD"])
        ename = "E"*string(Int64(v_["5"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "F1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(1.5))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["5"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "PBGBG"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(1.5))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["5"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "PBUBG"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(1.5))
        posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["6"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eI3PROD")
        arrset(ielftype,ie,iet_["eI3PROD"])
        ename = "E"*string(Int64(v_["6"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "F1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(1.5))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["6"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "F2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(1.5))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["6"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "PBGBG"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(1.5))
        posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["6"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "PBUBG"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(1.5))
        posev = findfirst(x->x=="V4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["7"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eI2PROD")
        arrset(ielftype,ie,iet_["eI2PROD"])
        ename = "E"*string(Int64(v_["7"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "F2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(1.5))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["7"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "PGBGB"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(1.5))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["7"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "PUBGB"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(1.5))
        posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["8"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eI3PROD")
        arrset(ielftype,ie,iet_["eI3PROD"])
        ename = "E"*string(Int64(v_["8"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "F1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(1.5))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["8"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "F2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(1.5))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["8"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "PGBGB"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(1.5))
        posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["8"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "PUBGB"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(1.5))
        posev = findfirst(x->x=="V4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["9"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en3PRODI")
        arrset(ielftype,ie,iet_["en3PRODI"])
        ename = "E"*string(Int64(v_["9"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "F1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(1.5))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["9"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "F2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(1.5))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["9"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "PBBBB"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(1.5))
        posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["9"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "PUBBB"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(1.5))
        posev = findfirst(x->x=="V4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E"*string(Int64(v_["9"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "PBUBB"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(1.5))
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
        loaset(pbm.grpar,ig,posgp,Float64(v_["RNBBUB"]))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               1138.416240
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        pb.clower[pb.nle+pb.neq+1:pb.m] = zeros(Float64,pb.nge)
        pb.cupper[1:pb.nge] = fill(Inf,pb.nge)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-COLR2-MY-9-4"
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
        f_   = IV_[1]*IV_[2]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = IV_[2]
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
        T = abs(GVAR_)
        SMALL = 1.0e-10
        LARGE = 1.0e+10
        ARG0 = T<=SMALL
        if ARG0
            FF = pbm.grpar[igr_][1]*log(SMALL)
        end
        if !ARG0
            FF = pbm.grpar[igr_][1]*log(T)
        end
        if ARG0
            GG = pbm.grpar[igr_][1]*LARGE
        end
        if !ARG0
            GG = pbm.grpar[igr_][1]/T
        end
        if ARG0
            HH = -pbm.grpar[igr_][1]*LARGE^2
        end
        if !ARG0
            HH = -pbm.grpar[igr_][1]/T^2
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

