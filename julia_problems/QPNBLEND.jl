function QPNBLEND(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : QPNBLEND
#    *********
# 
#    Source: a variant on the BLEND linear programming problem
#    with an additional NONCONVEX diagonal Hessian matrix as given by
#    N. I. M. Gould, "An algorithm for large-scale quadratic programming",
#    IMA J. Num. Anal (1991), 11, 299-324, problem class 4.
# 
#    SIF input: Nick Gould, January 1993
# 
#    classification = "C-QLR2-MN-83-74"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "QPNBLEND"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("1",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"1")
        ig,ig_,_ = s2mpj_ii("2",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"2")
        ig,ig_,_ = s2mpj_ii("3",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"3")
        ig,ig_,_ = s2mpj_ii("4",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"4")
        ig,ig_,_ = s2mpj_ii("5",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"5")
        ig,ig_,_ = s2mpj_ii("6",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"6")
        ig,ig_,_ = s2mpj_ii("7",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"7")
        ig,ig_,_ = s2mpj_ii("8",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"8")
        ig,ig_,_ = s2mpj_ii("9",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"9")
        ig,ig_,_ = s2mpj_ii("10",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"10")
        ig,ig_,_ = s2mpj_ii("11",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"11")
        ig,ig_,_ = s2mpj_ii("12",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"12")
        ig,ig_,_ = s2mpj_ii("13",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"13")
        ig,ig_,_ = s2mpj_ii("14",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"14")
        ig,ig_,_ = s2mpj_ii("15",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"15")
        ig,ig_,_ = s2mpj_ii("16",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"16")
        ig,ig_,_ = s2mpj_ii("17",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"17")
        ig,ig_,_ = s2mpj_ii("18",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"18")
        ig,ig_,_ = s2mpj_ii("19",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"19")
        ig,ig_,_ = s2mpj_ii("20",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"20")
        ig,ig_,_ = s2mpj_ii("21",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"21")
        ig,ig_,_ = s2mpj_ii("22",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"22")
        ig,ig_,_ = s2mpj_ii("23",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"23")
        ig,ig_,_ = s2mpj_ii("24",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"24")
        ig,ig_,_ = s2mpj_ii("25",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"25")
        ig,ig_,_ = s2mpj_ii("26",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"26")
        ig,ig_,_ = s2mpj_ii("27",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"27")
        ig,ig_,_ = s2mpj_ii("28",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"28")
        ig,ig_,_ = s2mpj_ii("29",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"29")
        ig,ig_,_ = s2mpj_ii("30",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"30")
        ig,ig_,_ = s2mpj_ii("31",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"31")
        ig,ig_,_ = s2mpj_ii("32",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"32")
        ig,ig_,_ = s2mpj_ii("33",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"33")
        ig,ig_,_ = s2mpj_ii("34",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"34")
        ig,ig_,_ = s2mpj_ii("35",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"35")
        ig,ig_,_ = s2mpj_ii("36",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"36")
        ig,ig_,_ = s2mpj_ii("37",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"37")
        ig,ig_,_ = s2mpj_ii("38",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"38")
        ig,ig_,_ = s2mpj_ii("39",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"39")
        ig,ig_,_ = s2mpj_ii("40",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"40")
        ig,ig_,_ = s2mpj_ii("41",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"41")
        ig,ig_,_ = s2mpj_ii("42",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"42")
        ig,ig_,_ = s2mpj_ii("43",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"43")
        ig,ig_,_ = s2mpj_ii("44",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"44")
        ig,ig_,_ = s2mpj_ii("45",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"45")
        ig,ig_,_ = s2mpj_ii("46",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"46")
        ig,ig_,_ = s2mpj_ii("47",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"47")
        ig,ig_,_ = s2mpj_ii("48",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"48")
        ig,ig_,_ = s2mpj_ii("49",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"49")
        ig,ig_,_ = s2mpj_ii("50",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"50")
        ig,ig_,_ = s2mpj_ii("51",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"51")
        ig,ig_,_ = s2mpj_ii("52",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"52")
        ig,ig_,_ = s2mpj_ii("53",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"53")
        ig,ig_,_ = s2mpj_ii("54",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"54")
        ig,ig_,_ = s2mpj_ii("55",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"55")
        ig,ig_,_ = s2mpj_ii("56",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"56")
        ig,ig_,_ = s2mpj_ii("57",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"57")
        ig,ig_,_ = s2mpj_ii("58",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"58")
        ig,ig_,_ = s2mpj_ii("59",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"59")
        ig,ig_,_ = s2mpj_ii("60",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"60")
        ig,ig_,_ = s2mpj_ii("61",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"61")
        ig,ig_,_ = s2mpj_ii("62",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"62")
        ig,ig_,_ = s2mpj_ii("63",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"63")
        ig,ig_,_ = s2mpj_ii("64",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"64")
        ig,ig_,_ = s2mpj_ii("65",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"65")
        ig,ig_,_ = s2mpj_ii("66",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"66")
        ig,ig_,_ = s2mpj_ii("67",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"67")
        ig,ig_,_ = s2mpj_ii("68",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"68")
        ig,ig_,_ = s2mpj_ii("69",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"69")
        ig,ig_,_ = s2mpj_ii("70",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"70")
        ig,ig_,_ = s2mpj_ii("71",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"71")
        ig,ig_,_ = s2mpj_ii("72",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"72")
        ig,ig_,_ = s2mpj_ii("73",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"73")
        ig,ig_,_ = s2mpj_ii("74",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"74")
        ig,ig_,_ = s2mpj_ii("C",ig_)
        arrset(gtype,ig,"<>")
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        ngrp   = length(ig_)
        iv,ix_,_ = s2mpj_ii("1",ix_)
        arrset(pb.xnames,iv,"1")
        ig = ig_["2"]
        pbm.A[ig,iv] += Float64(-.537)
        ig = ig_["3"]
        pbm.A[ig,iv] += Float64(-.131)
        iv,ix_,_ = s2mpj_ii("1",ix_)
        arrset(pb.xnames,iv,"1")
        ig = ig_["4"]
        pbm.A[ig,iv] += Float64(-.1155)
        ig = ig_["5"]
        pbm.A[ig,iv] += Float64(-.0365)
        iv,ix_,_ = s2mpj_ii("1",ix_)
        arrset(pb.xnames,iv,"1")
        ig = ig_["6"]
        pbm.A[ig,iv] += Float64(-.143)
        ig = ig_["7"]
        pbm.A[ig,iv] += Float64(-.037)
        iv,ix_,_ = s2mpj_ii("1",ix_)
        arrset(pb.xnames,iv,"1")
        ig = ig_["40"]
        pbm.A[ig,iv] += Float64(.003)
        ig = ig_["41"]
        pbm.A[ig,iv] += Float64(.0587)
        iv,ix_,_ = s2mpj_ii("1",ix_)
        arrset(pb.xnames,iv,"1")
        ig = ig_["42"]
        pbm.A[ig,iv] += Float64(.15)
        ig = ig_["43"]
        pbm.A[ig,iv] += Float64(.302)
        iv,ix_,_ = s2mpj_ii("1",ix_)
        arrset(pb.xnames,iv,"1")
        ig = ig_["67"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["C"]
        pbm.A[ig,iv] += Float64(3.2)
        iv,ix_,_ = s2mpj_ii("2",ix_)
        arrset(pb.xnames,iv,"2")
        ig = ig_["1"]
        pbm.A[ig,iv] += Float64(-.2931)
        ig = ig_["3"]
        pbm.A[ig,iv] += Float64(-.117)
        iv,ix_,_ = s2mpj_ii("2",ix_)
        arrset(pb.xnames,iv,"2")
        ig = ig_["4"]
        pbm.A[ig,iv] += Float64(-.0649)
        ig = ig_["5"]
        pbm.A[ig,iv] += Float64(-.1233)
        iv,ix_,_ = s2mpj_ii("2",ix_)
        arrset(pb.xnames,iv,"2")
        ig = ig_["6"]
        pbm.A[ig,iv] += Float64(-.2217)
        ig = ig_["8"]
        pbm.A[ig,iv] += Float64(-.18)
        iv,ix_,_ = s2mpj_ii("2",ix_)
        arrset(pb.xnames,iv,"2")
        ig = ig_["39"]
        pbm.A[ig,iv] += Float64(.0042)
        ig = ig_["40"]
        pbm.A[ig,iv] += Float64(.003)
        iv,ix_,_ = s2mpj_ii("2",ix_)
        arrset(pb.xnames,iv,"2")
        ig = ig_["41"]
        pbm.A[ig,iv] += Float64(.1053)
        ig = ig_["42"]
        pbm.A[ig,iv] += Float64(.185)
        iv,ix_,_ = s2mpj_ii("2",ix_)
        arrset(pb.xnames,iv,"2")
        ig = ig_["43"]
        pbm.A[ig,iv] += Float64(.384)
        ig = ig_["50"]
        pbm.A[ig,iv] += Float64(-.00862)
        iv,ix_,_ = s2mpj_ii("2",ix_)
        arrset(pb.xnames,iv,"2")
        ig = ig_["51"]
        pbm.A[ig,iv] += Float64(-.00862)
        ig = ig_["56"]
        pbm.A[ig,iv] += Float64(-.0101)
        iv,ix_,_ = s2mpj_ii("2",ix_)
        arrset(pb.xnames,iv,"2")
        ig = ig_["57"]
        pbm.A[ig,iv] += Float64(-.0101)
        ig = ig_["68"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("2",ix_)
        arrset(pb.xnames,iv,"2")
        ig = ig_["C"]
        pbm.A[ig,iv] += Float64(2.87)
        iv,ix_,_ = s2mpj_ii("3",ix_)
        arrset(pb.xnames,iv,"3")
        ig = ig_["2"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["9"]
        pbm.A[ig,iv] += Float64(-.0277)
        iv,ix_,_ = s2mpj_ii("3",ix_)
        arrset(pb.xnames,iv,"3")
        ig = ig_["10"]
        pbm.A[ig,iv] += Float64(-.0563)
        ig = ig_["11"]
        pbm.A[ig,iv] += Float64(-.199)
        iv,ix_,_ = s2mpj_ii("3",ix_)
        arrset(pb.xnames,iv,"3")
        ig = ig_["12"]
        pbm.A[ig,iv] += Float64(-.6873)
        ig = ig_["13"]
        pbm.A[ig,iv] += Float64(-.017)
        iv,ix_,_ = s2mpj_ii("3",ix_)
        arrset(pb.xnames,iv,"3")
        ig = ig_["40"]
        pbm.A[ig,iv] += Float64(.01303)
        ig = ig_["41"]
        pbm.A[ig,iv] += Float64(.0506)
        iv,ix_,_ = s2mpj_ii("3",ix_)
        arrset(pb.xnames,iv,"3")
        ig = ig_["42"]
        pbm.A[ig,iv] += Float64(.209)
        ig = ig_["43"]
        pbm.A[ig,iv] += Float64(.495)
        iv,ix_,_ = s2mpj_ii("3",ix_)
        arrset(pb.xnames,iv,"3")
        ig = ig_["65"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("4",ix_)
        arrset(pb.xnames,iv,"4")
        ig = ig_["1"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["9"]
        pbm.A[ig,iv] += Float64(-.0112)
        iv,ix_,_ = s2mpj_ii("4",ix_)
        arrset(pb.xnames,iv,"4")
        ig = ig_["10"]
        pbm.A[ig,iv] += Float64(-.0378)
        ig = ig_["11"]
        pbm.A[ig,iv] += Float64(-.1502)
        iv,ix_,_ = s2mpj_ii("4",ix_)
        arrset(pb.xnames,iv,"4")
        ig = ig_["12"]
        pbm.A[ig,iv] += Float64(-.7953)
        ig = ig_["13"]
        pbm.A[ig,iv] += Float64(-.0099)
        iv,ix_,_ = s2mpj_ii("4",ix_)
        arrset(pb.xnames,iv,"4")
        ig = ig_["40"]
        pbm.A[ig,iv] += Float64(.01303)
        ig = ig_["41"]
        pbm.A[ig,iv] += Float64(.0448)
        iv,ix_,_ = s2mpj_ii("4",ix_)
        arrset(pb.xnames,iv,"4")
        ig = ig_["42"]
        pbm.A[ig,iv] += Float64(.185)
        ig = ig_["43"]
        pbm.A[ig,iv] += Float64(.721)
        iv,ix_,_ = s2mpj_ii("4",ix_)
        arrset(pb.xnames,iv,"4")
        ig = ig_["65"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("5",ix_)
        arrset(pb.xnames,iv,"5")
        ig = ig_["9"]
        pbm.A[ig,iv] += Float64(-.175)
        ig = ig_["10"]
        pbm.A[ig,iv] += Float64(-.27)
        iv,ix_,_ = s2mpj_ii("5",ix_)
        arrset(pb.xnames,iv,"5")
        ig = ig_["11"]
        pbm.A[ig,iv] += Float64(-.028)
        ig = ig_["13"]
        pbm.A[ig,iv] += Float64(-.455)
        iv,ix_,_ = s2mpj_ii("5",ix_)
        arrset(pb.xnames,iv,"5")
        ig = ig_["21"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["40"]
        pbm.A[ig,iv] += Float64(.01303)
        iv,ix_,_ = s2mpj_ii("5",ix_)
        arrset(pb.xnames,iv,"5")
        ig = ig_["41"]
        pbm.A[ig,iv] += Float64(.0506)
        ig = ig_["42"]
        pbm.A[ig,iv] += Float64(.209)
        iv,ix_,_ = s2mpj_ii("5",ix_)
        arrset(pb.xnames,iv,"5")
        ig = ig_["43"]
        pbm.A[ig,iv] += Float64(.495)
        iv,ix_,_ = s2mpj_ii("6",ix_)
        arrset(pb.xnames,iv,"6")
        ig = ig_["9"]
        pbm.A[ig,iv] += Float64(-.271)
        ig = ig_["10"]
        pbm.A[ig,iv] += Float64(-.3285)
        iv,ix_,_ = s2mpj_ii("6",ix_)
        arrset(pb.xnames,iv,"6")
        ig = ig_["11"]
        pbm.A[ig,iv] += Float64(-.0255)
        ig = ig_["13"]
        pbm.A[ig,iv] += Float64(-.2656)
        iv,ix_,_ = s2mpj_ii("6",ix_)
        arrset(pb.xnames,iv,"6")
        ig = ig_["18"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["40"]
        pbm.A[ig,iv] += Float64(.01303)
        iv,ix_,_ = s2mpj_ii("6",ix_)
        arrset(pb.xnames,iv,"6")
        ig = ig_["41"]
        pbm.A[ig,iv] += Float64(.0506)
        ig = ig_["42"]
        pbm.A[ig,iv] += Float64(.209)
        iv,ix_,_ = s2mpj_ii("6",ix_)
        arrset(pb.xnames,iv,"6")
        ig = ig_["43"]
        pbm.A[ig,iv] += Float64(.495)
        iv,ix_,_ = s2mpj_ii("7",ix_)
        arrset(pb.xnames,iv,"7")
        ig = ig_["9"]
        pbm.A[ig,iv] += Float64(-.2836)
        ig = ig_["10"]
        pbm.A[ig,iv] += Float64(-.3285)
        iv,ix_,_ = s2mpj_ii("7",ix_)
        arrset(pb.xnames,iv,"7")
        ig = ig_["11"]
        pbm.A[ig,iv] += Float64(-.0241)
        ig = ig_["13"]
        pbm.A[ig,iv] += Float64(-.2502)
        iv,ix_,_ = s2mpj_ii("7",ix_)
        arrset(pb.xnames,iv,"7")
        ig = ig_["17"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["40"]
        pbm.A[ig,iv] += Float64(.01303)
        iv,ix_,_ = s2mpj_ii("7",ix_)
        arrset(pb.xnames,iv,"7")
        ig = ig_["41"]
        pbm.A[ig,iv] += Float64(.0506)
        ig = ig_["42"]
        pbm.A[ig,iv] += Float64(.209)
        iv,ix_,_ = s2mpj_ii("7",ix_)
        arrset(pb.xnames,iv,"7")
        ig = ig_["43"]
        pbm.A[ig,iv] += Float64(.495)
        iv,ix_,_ = s2mpj_ii("8",ix_)
        arrset(pb.xnames,iv,"8")
        ig = ig_["12"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["14"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("8",ix_)
        arrset(pb.xnames,iv,"8")
        ig = ig_["39"]
        pbm.A[ig,iv] += Float64(.0327)
        ig = ig_["41"]
        pbm.A[ig,iv] += Float64(.094)
        iv,ix_,_ = s2mpj_ii("8",ix_)
        arrset(pb.xnames,iv,"8")
        ig = ig_["42"]
        pbm.A[ig,iv] += Float64(.045)
        ig = ig_["43"]
        pbm.A[ig,iv] += Float64(.793)
        iv,ix_,_ = s2mpj_ii("8",ix_)
        arrset(pb.xnames,iv,"8")
        ig = ig_["C"]
        pbm.A[ig,iv] += Float64(.0044)
        iv,ix_,_ = s2mpj_ii("9",ix_)
        arrset(pb.xnames,iv,"9")
        ig = ig_["15"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["22"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("10",ix_)
        arrset(pb.xnames,iv,"10")
        ig = ig_["16"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["22"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("11",ix_)
        arrset(pb.xnames,iv,"11")
        ig = ig_["14"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["15"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("12",ix_)
        arrset(pb.xnames,iv,"12")
        ig = ig_["14"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["16"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("13",ix_)
        arrset(pb.xnames,iv,"13")
        ig = ig_["15"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["17"]
        pbm.A[ig,iv] += Float64(-.0588)
        iv,ix_,_ = s2mpj_ii("13",ix_)
        arrset(pb.xnames,iv,"13")
        ig = ig_["19"]
        pbm.A[ig,iv] += Float64(-.8145)
        ig = ig_["23"]
        pbm.A[ig,iv] += Float64(-.0091)
        iv,ix_,_ = s2mpj_ii("13",ix_)
        arrset(pb.xnames,iv,"13")
        ig = ig_["39"]
        pbm.A[ig,iv] += Float64(-.8239)
        ig = ig_["40"]
        pbm.A[ig,iv] += Float64(.0081)
        iv,ix_,_ = s2mpj_ii("13",ix_)
        arrset(pb.xnames,iv,"13")
        ig = ig_["41"]
        pbm.A[ig,iv] += Float64(-.2112)
        ig = ig_["42"]
        pbm.A[ig,iv] += Float64(.387)
        iv,ix_,_ = s2mpj_ii("13",ix_)
        arrset(pb.xnames,iv,"13")
        ig = ig_["43"]
        pbm.A[ig,iv] += Float64(1.03)
        ig = ig_["69"]
        pbm.A[ig,iv] += Float64(1.3)
        iv,ix_,_ = s2mpj_ii("13",ix_)
        arrset(pb.xnames,iv,"13")
        ig = ig_["C"]
        pbm.A[ig,iv] += Float64(.07)
        iv,ix_,_ = s2mpj_ii("14",ix_)
        arrset(pb.xnames,iv,"14")
        ig = ig_["16"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["18"]
        pbm.A[ig,iv] += Float64(-.0404)
        iv,ix_,_ = s2mpj_ii("14",ix_)
        arrset(pb.xnames,iv,"14")
        ig = ig_["20"]
        pbm.A[ig,iv] += Float64(-.8564)
        ig = ig_["23"]
        pbm.A[ig,iv] += Float64(-.0069)
        iv,ix_,_ = s2mpj_ii("14",ix_)
        arrset(pb.xnames,iv,"14")
        ig = ig_["39"]
        pbm.A[ig,iv] += Float64(-.7689)
        ig = ig_["40"]
        pbm.A[ig,iv] += Float64(.0063)
        iv,ix_,_ = s2mpj_ii("14",ix_)
        arrset(pb.xnames,iv,"14")
        ig = ig_["41"]
        pbm.A[ig,iv] += Float64(-.156)
        ig = ig_["42"]
        pbm.A[ig,iv] += Float64(.297)
        iv,ix_,_ = s2mpj_ii("14",ix_)
        arrset(pb.xnames,iv,"14")
        ig = ig_["43"]
        pbm.A[ig,iv] += Float64(.792)
        ig = ig_["69"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("14",ix_)
        arrset(pb.xnames,iv,"14")
        ig = ig_["C"]
        pbm.A[ig,iv] += Float64(.0378)
        iv,ix_,_ = s2mpj_ii("15",ix_)
        arrset(pb.xnames,iv,"15")
        ig = ig_["5"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["21"]
        pbm.A[ig,iv] += Float64(-.3321)
        iv,ix_,_ = s2mpj_ii("15",ix_)
        arrset(pb.xnames,iv,"15")
        ig = ig_["22"]
        pbm.A[ig,iv] += Float64(-.5875)
        ig = ig_["23"]
        pbm.A[ig,iv] += Float64(-.362)
        iv,ix_,_ = s2mpj_ii("15",ix_)
        arrset(pb.xnames,iv,"15")
        ig = ig_["39"]
        pbm.A[ig,iv] += Float64(2.3)
        ig = ig_["41"]
        pbm.A[ig,iv] += Float64(-.2049)
        iv,ix_,_ = s2mpj_ii("15",ix_)
        arrset(pb.xnames,iv,"15")
        ig = ig_["42"]
        pbm.A[ig,iv] += Float64(.826)
        ig = ig_["43"]
        pbm.A[ig,iv] += Float64(14.61)
        iv,ix_,_ = s2mpj_ii("15",ix_)
        arrset(pb.xnames,iv,"15")
        ig = ig_["65"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["70"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("15",ix_)
        arrset(pb.xnames,iv,"15")
        ig = ig_["C"]
        pbm.A[ig,iv] += Float64(.155)
        iv,ix_,_ = s2mpj_ii("16",ix_)
        arrset(pb.xnames,iv,"16")
        ig = ig_["6"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["21"]
        pbm.A[ig,iv] += Float64(-.3321)
        iv,ix_,_ = s2mpj_ii("16",ix_)
        arrset(pb.xnames,iv,"16")
        ig = ig_["22"]
        pbm.A[ig,iv] += Float64(-.5875)
        ig = ig_["23"]
        pbm.A[ig,iv] += Float64(-.362)
        iv,ix_,_ = s2mpj_ii("16",ix_)
        arrset(pb.xnames,iv,"16")
        ig = ig_["39"]
        pbm.A[ig,iv] += Float64(2.3)
        ig = ig_["41"]
        pbm.A[ig,iv] += Float64(-.2049)
        iv,ix_,_ = s2mpj_ii("16",ix_)
        arrset(pb.xnames,iv,"16")
        ig = ig_["42"]
        pbm.A[ig,iv] += Float64(.826)
        ig = ig_["43"]
        pbm.A[ig,iv] += Float64(14.61)
        iv,ix_,_ = s2mpj_ii("16",ix_)
        arrset(pb.xnames,iv,"16")
        ig = ig_["66"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["70"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("16",ix_)
        arrset(pb.xnames,iv,"16")
        ig = ig_["C"]
        pbm.A[ig,iv] += Float64(.155)
        iv,ix_,_ = s2mpj_ii("17",ix_)
        arrset(pb.xnames,iv,"17")
        ig = ig_["4"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["21"]
        pbm.A[ig,iv] += Float64(-.2414)
        iv,ix_,_ = s2mpj_ii("17",ix_)
        arrset(pb.xnames,iv,"17")
        ig = ig_["22"]
        pbm.A[ig,iv] += Float64(-.6627)
        ig = ig_["23"]
        pbm.A[ig,iv] += Float64(-.293)
        iv,ix_,_ = s2mpj_ii("17",ix_)
        arrset(pb.xnames,iv,"17")
        ig = ig_["39"]
        pbm.A[ig,iv] += Float64(2.3)
        ig = ig_["41"]
        pbm.A[ig,iv] += Float64(-.1531)
        iv,ix_,_ = s2mpj_ii("17",ix_)
        arrset(pb.xnames,iv,"17")
        ig = ig_["42"]
        pbm.A[ig,iv] += Float64(.826)
        ig = ig_["43"]
        pbm.A[ig,iv] += Float64(14.61)
        iv,ix_,_ = s2mpj_ii("17",ix_)
        arrset(pb.xnames,iv,"17")
        ig = ig_["65"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["70"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("17",ix_)
        arrset(pb.xnames,iv,"17")
        ig = ig_["C"]
        pbm.A[ig,iv] += Float64(.155)
        iv,ix_,_ = s2mpj_ii("18",ix_)
        arrset(pb.xnames,iv,"18")
        ig = ig_["21"]
        pbm.A[ig,iv] += Float64(-.2414)
        ig = ig_["22"]
        pbm.A[ig,iv] += Float64(-.6627)
        iv,ix_,_ = s2mpj_ii("18",ix_)
        arrset(pb.xnames,iv,"18")
        ig = ig_["23"]
        pbm.A[ig,iv] += Float64(-.293)
        ig = ig_["28"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("18",ix_)
        arrset(pb.xnames,iv,"18")
        ig = ig_["39"]
        pbm.A[ig,iv] += Float64(2.3)
        ig = ig_["41"]
        pbm.A[ig,iv] += Float64(-.1531)
        iv,ix_,_ = s2mpj_ii("18",ix_)
        arrset(pb.xnames,iv,"18")
        ig = ig_["42"]
        pbm.A[ig,iv] += Float64(.826)
        ig = ig_["43"]
        pbm.A[ig,iv] += Float64(14.61)
        iv,ix_,_ = s2mpj_ii("18",ix_)
        arrset(pb.xnames,iv,"18")
        ig = ig_["70"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["C"]
        pbm.A[ig,iv] += Float64(.155)
        iv,ix_,_ = s2mpj_ii("19",ix_)
        arrset(pb.xnames,iv,"19")
        ig = ig_["5"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["10"]
        pbm.A[ig,iv] += Float64(-.0185)
        iv,ix_,_ = s2mpj_ii("19",ix_)
        arrset(pb.xnames,iv,"19")
        ig = ig_["13"]
        pbm.A[ig,iv] += Float64(-.0568)
        ig = ig_["24"]
        pbm.A[ig,iv] += Float64(-.0806)
        iv,ix_,_ = s2mpj_ii("19",ix_)
        arrset(pb.xnames,iv,"19")
        ig = ig_["25"]
        pbm.A[ig,iv] += Float64(-.0658)
        ig = ig_["26"]
        pbm.A[ig,iv] += Float64(-.0328)
        iv,ix_,_ = s2mpj_ii("19",ix_)
        arrset(pb.xnames,iv,"19")
        ig = ig_["27"]
        pbm.A[ig,iv] += Float64(-.4934)
        ig = ig_["28"]
        pbm.A[ig,iv] += Float64(-.2922)
        iv,ix_,_ = s2mpj_ii("19",ix_)
        arrset(pb.xnames,iv,"19")
        ig = ig_["29"]
        pbm.A[ig,iv] += Float64(-.0096)
        ig = ig_["40"]
        pbm.A[ig,iv] += Float64(-.0654)
        iv,ix_,_ = s2mpj_ii("19",ix_)
        arrset(pb.xnames,iv,"19")
        ig = ig_["41"]
        pbm.A[ig,iv] += Float64(-.2535)
        ig = ig_["42"]
        pbm.A[ig,iv] += Float64(.632)
        iv,ix_,_ = s2mpj_ii("19",ix_)
        arrset(pb.xnames,iv,"19")
        ig = ig_["43"]
        pbm.A[ig,iv] += Float64(.6807)
        ig = ig_["65"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("19",ix_)
        arrset(pb.xnames,iv,"19")
        ig = ig_["71"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["C"]
        pbm.A[ig,iv] += Float64(.0528)
        iv,ix_,_ = s2mpj_ii("20",ix_)
        arrset(pb.xnames,iv,"20")
        ig = ig_["6"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["10"]
        pbm.A[ig,iv] += Float64(-.0185)
        iv,ix_,_ = s2mpj_ii("20",ix_)
        arrset(pb.xnames,iv,"20")
        ig = ig_["13"]
        pbm.A[ig,iv] += Float64(-.0568)
        ig = ig_["24"]
        pbm.A[ig,iv] += Float64(-.0806)
        iv,ix_,_ = s2mpj_ii("20",ix_)
        arrset(pb.xnames,iv,"20")
        ig = ig_["25"]
        pbm.A[ig,iv] += Float64(-.0658)
        ig = ig_["26"]
        pbm.A[ig,iv] += Float64(-.0328)
        iv,ix_,_ = s2mpj_ii("20",ix_)
        arrset(pb.xnames,iv,"20")
        ig = ig_["27"]
        pbm.A[ig,iv] += Float64(-.4934)
        ig = ig_["28"]
        pbm.A[ig,iv] += Float64(-.2922)
        iv,ix_,_ = s2mpj_ii("20",ix_)
        arrset(pb.xnames,iv,"20")
        ig = ig_["29"]
        pbm.A[ig,iv] += Float64(-.0096)
        ig = ig_["40"]
        pbm.A[ig,iv] += Float64(-.0654)
        iv,ix_,_ = s2mpj_ii("20",ix_)
        arrset(pb.xnames,iv,"20")
        ig = ig_["41"]
        pbm.A[ig,iv] += Float64(-.2535)
        ig = ig_["42"]
        pbm.A[ig,iv] += Float64(.632)
        iv,ix_,_ = s2mpj_ii("20",ix_)
        arrset(pb.xnames,iv,"20")
        ig = ig_["43"]
        pbm.A[ig,iv] += Float64(.6807)
        ig = ig_["66"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("20",ix_)
        arrset(pb.xnames,iv,"20")
        ig = ig_["71"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["C"]
        pbm.A[ig,iv] += Float64(.0528)
        iv,ix_,_ = s2mpj_ii("21",ix_)
        arrset(pb.xnames,iv,"21")
        ig = ig_["4"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["10"]
        pbm.A[ig,iv] += Float64(-.0184)
        iv,ix_,_ = s2mpj_ii("21",ix_)
        arrset(pb.xnames,iv,"21")
        ig = ig_["13"]
        pbm.A[ig,iv] += Float64(-.0564)
        ig = ig_["24"]
        pbm.A[ig,iv] += Float64(-.078)
        iv,ix_,_ = s2mpj_ii("21",ix_)
        arrset(pb.xnames,iv,"21")
        ig = ig_["25"]
        pbm.A[ig,iv] += Float64(-.0655)
        ig = ig_["26"]
        pbm.A[ig,iv] += Float64(-.0303)
        iv,ix_,_ = s2mpj_ii("21",ix_)
        arrset(pb.xnames,iv,"21")
        ig = ig_["27"]
        pbm.A[ig,iv] += Float64(-.475)
        ig = ig_["28"]
        pbm.A[ig,iv] += Float64(-.305)
        iv,ix_,_ = s2mpj_ii("21",ix_)
        arrset(pb.xnames,iv,"21")
        ig = ig_["40"]
        pbm.A[ig,iv] += Float64(-.0654)
        ig = ig_["41"]
        pbm.A[ig,iv] += Float64(-.2703)
        iv,ix_,_ = s2mpj_ii("21",ix_)
        arrset(pb.xnames,iv,"21")
        ig = ig_["42"]
        pbm.A[ig,iv] += Float64(.632)
        ig = ig_["43"]
        pbm.A[ig,iv] += Float64(.6807)
        iv,ix_,_ = s2mpj_ii("21",ix_)
        arrset(pb.xnames,iv,"21")
        ig = ig_["65"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["71"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("21",ix_)
        arrset(pb.xnames,iv,"21")
        ig = ig_["C"]
        pbm.A[ig,iv] += Float64(.0528)
        iv,ix_,_ = s2mpj_ii("22",ix_)
        arrset(pb.xnames,iv,"22")
        ig = ig_["3"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["10"]
        pbm.A[ig,iv] += Float64(-.0184)
        iv,ix_,_ = s2mpj_ii("22",ix_)
        arrset(pb.xnames,iv,"22")
        ig = ig_["13"]
        pbm.A[ig,iv] += Float64(-.0564)
        ig = ig_["24"]
        pbm.A[ig,iv] += Float64(-.078)
        iv,ix_,_ = s2mpj_ii("22",ix_)
        arrset(pb.xnames,iv,"22")
        ig = ig_["25"]
        pbm.A[ig,iv] += Float64(-.0655)
        ig = ig_["26"]
        pbm.A[ig,iv] += Float64(-.0303)
        iv,ix_,_ = s2mpj_ii("22",ix_)
        arrset(pb.xnames,iv,"22")
        ig = ig_["27"]
        pbm.A[ig,iv] += Float64(-.475)
        ig = ig_["28"]
        pbm.A[ig,iv] += Float64(-.305)
        iv,ix_,_ = s2mpj_ii("22",ix_)
        arrset(pb.xnames,iv,"22")
        ig = ig_["40"]
        pbm.A[ig,iv] += Float64(-.0654)
        ig = ig_["41"]
        pbm.A[ig,iv] += Float64(-.2703)
        iv,ix_,_ = s2mpj_ii("22",ix_)
        arrset(pb.xnames,iv,"22")
        ig = ig_["42"]
        pbm.A[ig,iv] += Float64(.632)
        ig = ig_["43"]
        pbm.A[ig,iv] += Float64(.6807)
        iv,ix_,_ = s2mpj_ii("22",ix_)
        arrset(pb.xnames,iv,"22")
        ig = ig_["65"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["71"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("22",ix_)
        arrset(pb.xnames,iv,"22")
        ig = ig_["C"]
        pbm.A[ig,iv] += Float64(.0528)
        iv,ix_,_ = s2mpj_ii("23",ix_)
        arrset(pb.xnames,iv,"23")
        ig = ig_["13"]
        pbm.A[ig,iv] += Float64(.76)
        ig = ig_["25"]
        pbm.A[ig,iv] += Float64(.5714)
        iv,ix_,_ = s2mpj_ii("23",ix_)
        arrset(pb.xnames,iv,"23")
        ig = ig_["30"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["40"]
        pbm.A[ig,iv] += Float64(.1869)
        iv,ix_,_ = s2mpj_ii("23",ix_)
        arrset(pb.xnames,iv,"23")
        ig = ig_["41"]
        pbm.A[ig,iv] += Float64(.2796)
        ig = ig_["42"]
        pbm.A[ig,iv] += Float64(2.241)
        iv,ix_,_ = s2mpj_ii("23",ix_)
        arrset(pb.xnames,iv,"23")
        ig = ig_["43"]
        pbm.A[ig,iv] += Float64(2.766)
        ig = ig_["72"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("23",ix_)
        arrset(pb.xnames,iv,"23")
        ig = ig_["C"]
        pbm.A[ig,iv] += Float64(.128)
        iv,ix_,_ = s2mpj_ii("24",ix_)
        arrset(pb.xnames,iv,"24")
        ig = ig_["9"]
        pbm.A[ig,iv] += Float64(-.0571)
        ig = ig_["10"]
        pbm.A[ig,iv] += Float64(-.0114)
        iv,ix_,_ = s2mpj_ii("24",ix_)
        arrset(pb.xnames,iv,"24")
        ig = ig_["13"]
        pbm.A[ig,iv] += Float64(.6571)
        ig = ig_["24"]
        pbm.A[ig,iv] += Float64(.5714)
        iv,ix_,_ = s2mpj_ii("24",ix_)
        arrset(pb.xnames,iv,"24")
        ig = ig_["31"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["40"]
        pbm.A[ig,iv] += Float64(.1724)
        iv,ix_,_ = s2mpj_ii("24",ix_)
        arrset(pb.xnames,iv,"24")
        ig = ig_["41"]
        pbm.A[ig,iv] += Float64(.2579)
        ig = ig_["42"]
        pbm.A[ig,iv] += Float64(2.067)
        iv,ix_,_ = s2mpj_ii("24",ix_)
        arrset(pb.xnames,iv,"24")
        ig = ig_["43"]
        pbm.A[ig,iv] += Float64(2.552)
        ig = ig_["72"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("24",ix_)
        arrset(pb.xnames,iv,"24")
        ig = ig_["C"]
        pbm.A[ig,iv] += Float64(.118)
        iv,ix_,_ = s2mpj_ii("25",ix_)
        arrset(pb.xnames,iv,"25")
        ig = ig_["9"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["25"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("26",ix_)
        arrset(pb.xnames,iv,"26")
        ig = ig_["10"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["24"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("27",ix_)
        arrset(pb.xnames,iv,"27")
        ig = ig_["10"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["13"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("28",ix_)
        arrset(pb.xnames,iv,"28")
        ig = ig_["11"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["32"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("28",ix_)
        arrset(pb.xnames,iv,"28")
        ig = ig_["44"]
        pbm.A[ig,iv] += Float64(-7.95)
        ig = ig_["45"]
        pbm.A[ig,iv] += Float64(-8.7)
        iv,ix_,_ = s2mpj_ii("28",ix_)
        arrset(pb.xnames,iv,"28")
        ig = ig_["46"]
        pbm.A[ig,iv] += Float64(-3.0)
        ig = ig_["47"]
        pbm.A[ig,iv] += Float64(14.0)
        iv,ix_,_ = s2mpj_ii("28",ix_)
        arrset(pb.xnames,iv,"28")
        ig = ig_["48"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["49"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("29",ix_)
        arrset(pb.xnames,iv,"29")
        ig = ig_["23"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["32"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("29",ix_)
        arrset(pb.xnames,iv,"29")
        ig = ig_["44"]
        pbm.A[ig,iv] += Float64(-8.84)
        ig = ig_["45"]
        pbm.A[ig,iv] += Float64(-9.45)
        iv,ix_,_ = s2mpj_ii("29",ix_)
        arrset(pb.xnames,iv,"29")
        ig = ig_["46"]
        pbm.A[ig,iv] += Float64(-3.0)
        ig = ig_["47"]
        pbm.A[ig,iv] += Float64(12.0)
        iv,ix_,_ = s2mpj_ii("29",ix_)
        arrset(pb.xnames,iv,"29")
        ig = ig_["48"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["49"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("30",ix_)
        arrset(pb.xnames,iv,"30")
        ig = ig_["19"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["32"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("30",ix_)
        arrset(pb.xnames,iv,"30")
        ig = ig_["44"]
        pbm.A[ig,iv] += Float64(-9.43)
        ig = ig_["45"]
        pbm.A[ig,iv] += Float64(-9.57)
        iv,ix_,_ = s2mpj_ii("30",ix_)
        arrset(pb.xnames,iv,"30")
        ig = ig_["46"]
        pbm.A[ig,iv] += Float64(-3.0)
        ig = ig_["47"]
        pbm.A[ig,iv] += Float64(3.5)
        iv,ix_,_ = s2mpj_ii("30",ix_)
        arrset(pb.xnames,iv,"30")
        ig = ig_["48"]
        pbm.A[ig,iv] += Float64(.233)
        ig = ig_["49"]
        pbm.A[ig,iv] += Float64(-.358)
        iv,ix_,_ = s2mpj_ii("31",ix_)
        arrset(pb.xnames,iv,"31")
        ig = ig_["20"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["32"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("31",ix_)
        arrset(pb.xnames,iv,"31")
        ig = ig_["44"]
        pbm.A[ig,iv] += Float64(-9.03)
        ig = ig_["45"]
        pbm.A[ig,iv] += Float64(-9.32)
        iv,ix_,_ = s2mpj_ii("31",ix_)
        arrset(pb.xnames,iv,"31")
        ig = ig_["46"]
        pbm.A[ig,iv] += Float64(-3.0)
        ig = ig_["47"]
        pbm.A[ig,iv] += Float64(3.5)
        iv,ix_,_ = s2mpj_ii("31",ix_)
        arrset(pb.xnames,iv,"31")
        ig = ig_["48"]
        pbm.A[ig,iv] += Float64(.205)
        ig = ig_["49"]
        pbm.A[ig,iv] += Float64(-.333)
        iv,ix_,_ = s2mpj_ii("32",ix_)
        arrset(pb.xnames,iv,"32")
        ig = ig_["27"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["32"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("32",ix_)
        arrset(pb.xnames,iv,"32")
        ig = ig_["44"]
        pbm.A[ig,iv] += Float64(-9.23)
        ig = ig_["45"]
        pbm.A[ig,iv] += Float64(-9.22)
        iv,ix_,_ = s2mpj_ii("32",ix_)
        arrset(pb.xnames,iv,"32")
        ig = ig_["46"]
        pbm.A[ig,iv] += Float64(-3.0)
        ig = ig_["47"]
        pbm.A[ig,iv] += Float64(6.0)
        iv,ix_,_ = s2mpj_ii("32",ix_)
        arrset(pb.xnames,iv,"32")
        ig = ig_["48"]
        pbm.A[ig,iv] += Float64(.381)
        ig = ig_["49"]
        pbm.A[ig,iv] += Float64(-.509)
        iv,ix_,_ = s2mpj_ii("33",ix_)
        arrset(pb.xnames,iv,"33")
        ig = ig_["30"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["32"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("33",ix_)
        arrset(pb.xnames,iv,"33")
        ig = ig_["44"]
        pbm.A[ig,iv] += Float64(-9.4)
        ig = ig_["45"]
        pbm.A[ig,iv] += Float64(-9.85)
        iv,ix_,_ = s2mpj_ii("33",ix_)
        arrset(pb.xnames,iv,"33")
        ig = ig_["46"]
        pbm.A[ig,iv] += Float64(-3.0)
        ig = ig_["47"]
        pbm.A[ig,iv] += Float64(2.5)
        iv,ix_,_ = s2mpj_ii("33",ix_)
        arrset(pb.xnames,iv,"33")
        ig = ig_["48"]
        pbm.A[ig,iv] += Float64(.39)
        ig = ig_["49"]
        pbm.A[ig,iv] += Float64(-.77)
        iv,ix_,_ = s2mpj_ii("34",ix_)
        arrset(pb.xnames,iv,"34")
        ig = ig_["31"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["32"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("34",ix_)
        arrset(pb.xnames,iv,"34")
        ig = ig_["44"]
        pbm.A[ig,iv] += Float64(-9.74)
        ig = ig_["45"]
        pbm.A[ig,iv] += Float64(-10.1)
        iv,ix_,_ = s2mpj_ii("34",ix_)
        arrset(pb.xnames,iv,"34")
        ig = ig_["46"]
        pbm.A[ig,iv] += Float64(-3.0)
        ig = ig_["47"]
        pbm.A[ig,iv] += Float64(3.3)
        iv,ix_,_ = s2mpj_ii("34",ix_)
        arrset(pb.xnames,iv,"34")
        ig = ig_["48"]
        pbm.A[ig,iv] += Float64(.233)
        ig = ig_["49"]
        pbm.A[ig,iv] += Float64(-.58)
        iv,ix_,_ = s2mpj_ii("35",ix_)
        arrset(pb.xnames,iv,"35")
        ig = ig_["10"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["32"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("35",ix_)
        arrset(pb.xnames,iv,"35")
        ig = ig_["44"]
        pbm.A[ig,iv] += Float64(-9.74)
        ig = ig_["45"]
        pbm.A[ig,iv] += Float64(-9.9)
        iv,ix_,_ = s2mpj_ii("35",ix_)
        arrset(pb.xnames,iv,"35")
        ig = ig_["46"]
        pbm.A[ig,iv] += Float64(-3.0)
        ig = ig_["47"]
        pbm.A[ig,iv] += Float64(66.0)
        iv,ix_,_ = s2mpj_ii("35",ix_)
        arrset(pb.xnames,iv,"35")
        ig = ig_["48"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["49"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("36",ix_)
        arrset(pb.xnames,iv,"36")
        ig = ig_["44"]
        pbm.A[ig,iv] += Float64(-.493)
        ig = ig_["45"]
        pbm.A[ig,iv] += Float64(-.165)
        iv,ix_,_ = s2mpj_ii("36",ix_)
        arrset(pb.xnames,iv,"36")
        ig = ig_["46"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["C"]
        pbm.A[ig,iv] += Float64(.0924)
        iv,ix_,_ = s2mpj_ii("37",ix_)
        arrset(pb.xnames,iv,"37")
        ig = ig_["32"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["44"]
        pbm.A[ig,iv] += Float64(10.03)
        iv,ix_,_ = s2mpj_ii("37",ix_)
        arrset(pb.xnames,iv,"37")
        ig = ig_["45"]
        pbm.A[ig,iv] += Float64(10.03)
        ig = ig_["47"]
        pbm.A[ig,iv] += Float64(-9.5)
        iv,ix_,_ = s2mpj_ii("37",ix_)
        arrset(pb.xnames,iv,"37")
        ig = ig_["48"]
        pbm.A[ig,iv] += Float64(-.5)
        ig = ig_["49"]
        pbm.A[ig,iv] += Float64(.5)
        iv,ix_,_ = s2mpj_ii("37",ix_)
        arrset(pb.xnames,iv,"37")
        ig = ig_["73"]
        pbm.A[ig,iv] += Float64(.64)
        ig = ig_["74"]
        pbm.A[ig,iv] += Float64(.35)
        iv,ix_,_ = s2mpj_ii("37",ix_)
        arrset(pb.xnames,iv,"37")
        ig = ig_["C"]
        pbm.A[ig,iv] += Float64(-5.36)
        iv,ix_,_ = s2mpj_ii("38",ix_)
        arrset(pb.xnames,iv,"38")
        ig = ig_["11"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["33"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("38",ix_)
        arrset(pb.xnames,iv,"38")
        ig = ig_["50"]
        pbm.A[ig,iv] += Float64(-7.98)
        ig = ig_["51"]
        pbm.A[ig,iv] += Float64(-8.58)
        iv,ix_,_ = s2mpj_ii("38",ix_)
        arrset(pb.xnames,iv,"38")
        ig = ig_["52"]
        pbm.A[ig,iv] += Float64(-3.0)
        ig = ig_["53"]
        pbm.A[ig,iv] += Float64(14.0)
        iv,ix_,_ = s2mpj_ii("38",ix_)
        arrset(pb.xnames,iv,"38")
        ig = ig_["54"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["55"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("39",ix_)
        arrset(pb.xnames,iv,"39")
        ig = ig_["23"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["33"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("39",ix_)
        arrset(pb.xnames,iv,"39")
        ig = ig_["50"]
        pbm.A[ig,iv] += Float64(-8.87)
        ig = ig_["51"]
        pbm.A[ig,iv] += Float64(-9.33)
        iv,ix_,_ = s2mpj_ii("39",ix_)
        arrset(pb.xnames,iv,"39")
        ig = ig_["52"]
        pbm.A[ig,iv] += Float64(-3.0)
        ig = ig_["53"]
        pbm.A[ig,iv] += Float64(12.0)
        iv,ix_,_ = s2mpj_ii("39",ix_)
        arrset(pb.xnames,iv,"39")
        ig = ig_["54"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["55"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("40",ix_)
        arrset(pb.xnames,iv,"40")
        ig = ig_["19"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["33"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("40",ix_)
        arrset(pb.xnames,iv,"40")
        ig = ig_["50"]
        pbm.A[ig,iv] += Float64(-9.46)
        ig = ig_["51"]
        pbm.A[ig,iv] += Float64(-9.45)
        iv,ix_,_ = s2mpj_ii("40",ix_)
        arrset(pb.xnames,iv,"40")
        ig = ig_["52"]
        pbm.A[ig,iv] += Float64(-3.0)
        ig = ig_["53"]
        pbm.A[ig,iv] += Float64(3.5)
        iv,ix_,_ = s2mpj_ii("40",ix_)
        arrset(pb.xnames,iv,"40")
        ig = ig_["54"]
        pbm.A[ig,iv] += Float64(.233)
        ig = ig_["55"]
        pbm.A[ig,iv] += Float64(-.358)
        iv,ix_,_ = s2mpj_ii("41",ix_)
        arrset(pb.xnames,iv,"41")
        ig = ig_["20"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["33"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("41",ix_)
        arrset(pb.xnames,iv,"41")
        ig = ig_["50"]
        pbm.A[ig,iv] += Float64(-9.06)
        ig = ig_["51"]
        pbm.A[ig,iv] += Float64(-9.2)
        iv,ix_,_ = s2mpj_ii("41",ix_)
        arrset(pb.xnames,iv,"41")
        ig = ig_["52"]
        pbm.A[ig,iv] += Float64(-3.0)
        ig = ig_["53"]
        pbm.A[ig,iv] += Float64(3.5)
        iv,ix_,_ = s2mpj_ii("41",ix_)
        arrset(pb.xnames,iv,"41")
        ig = ig_["54"]
        pbm.A[ig,iv] += Float64(.205)
        ig = ig_["55"]
        pbm.A[ig,iv] += Float64(-.333)
        iv,ix_,_ = s2mpj_ii("42",ix_)
        arrset(pb.xnames,iv,"42")
        ig = ig_["27"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["33"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("42",ix_)
        arrset(pb.xnames,iv,"42")
        ig = ig_["50"]
        pbm.A[ig,iv] += Float64(-9.26)
        ig = ig_["51"]
        pbm.A[ig,iv] += Float64(-9.13)
        iv,ix_,_ = s2mpj_ii("42",ix_)
        arrset(pb.xnames,iv,"42")
        ig = ig_["52"]
        pbm.A[ig,iv] += Float64(-3.0)
        ig = ig_["53"]
        pbm.A[ig,iv] += Float64(6.0)
        iv,ix_,_ = s2mpj_ii("42",ix_)
        arrset(pb.xnames,iv,"42")
        ig = ig_["54"]
        pbm.A[ig,iv] += Float64(.318)
        ig = ig_["55"]
        pbm.A[ig,iv] += Float64(-.509)
        iv,ix_,_ = s2mpj_ii("43",ix_)
        arrset(pb.xnames,iv,"43")
        ig = ig_["10"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["33"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("43",ix_)
        arrset(pb.xnames,iv,"43")
        ig = ig_["50"]
        pbm.A[ig,iv] += Float64(-9.77)
        ig = ig_["51"]
        pbm.A[ig,iv] += Float64(-9.78)
        iv,ix_,_ = s2mpj_ii("43",ix_)
        arrset(pb.xnames,iv,"43")
        ig = ig_["52"]
        pbm.A[ig,iv] += Float64(-3.0)
        ig = ig_["53"]
        pbm.A[ig,iv] += Float64(66.0)
        iv,ix_,_ = s2mpj_ii("43",ix_)
        arrset(pb.xnames,iv,"43")
        ig = ig_["54"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["55"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("44",ix_)
        arrset(pb.xnames,iv,"44")
        ig = ig_["50"]
        pbm.A[ig,iv] += Float64(-.435)
        ig = ig_["51"]
        pbm.A[ig,iv] += Float64(-.208)
        iv,ix_,_ = s2mpj_ii("44",ix_)
        arrset(pb.xnames,iv,"44")
        ig = ig_["52"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["C"]
        pbm.A[ig,iv] += Float64(.0924)
        iv,ix_,_ = s2mpj_ii("45",ix_)
        arrset(pb.xnames,iv,"45")
        ig = ig_["33"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["50"]
        pbm.A[ig,iv] += Float64(9.65)
        iv,ix_,_ = s2mpj_ii("45",ix_)
        arrset(pb.xnames,iv,"45")
        ig = ig_["51"]
        pbm.A[ig,iv] += Float64(9.65)
        ig = ig_["53"]
        pbm.A[ig,iv] += Float64(-9.5)
        iv,ix_,_ = s2mpj_ii("45",ix_)
        arrset(pb.xnames,iv,"45")
        ig = ig_["54"]
        pbm.A[ig,iv] += Float64(-.5)
        ig = ig_["55"]
        pbm.A[ig,iv] += Float64(.5)
        iv,ix_,_ = s2mpj_ii("45",ix_)
        arrset(pb.xnames,iv,"45")
        ig = ig_["73"]
        pbm.A[ig,iv] += Float64(-.36)
        ig = ig_["74"]
        pbm.A[ig,iv] += Float64(.35)
        iv,ix_,_ = s2mpj_ii("45",ix_)
        arrset(pb.xnames,iv,"45")
        ig = ig_["C"]
        pbm.A[ig,iv] += Float64(-5.08)
        iv,ix_,_ = s2mpj_ii("46",ix_)
        arrset(pb.xnames,iv,"46")
        ig = ig_["11"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["36"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("46",ix_)
        arrset(pb.xnames,iv,"46")
        ig = ig_["56"]
        pbm.A[ig,iv] += Float64(-7.99)
        ig = ig_["57"]
        pbm.A[ig,iv] += Float64(-8.59)
        iv,ix_,_ = s2mpj_ii("46",ix_)
        arrset(pb.xnames,iv,"46")
        ig = ig_["58"]
        pbm.A[ig,iv] += Float64(-3.0)
        ig = ig_["59"]
        pbm.A[ig,iv] += Float64(14.0)
        iv,ix_,_ = s2mpj_ii("46",ix_)
        arrset(pb.xnames,iv,"46")
        ig = ig_["60"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["61"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("47",ix_)
        arrset(pb.xnames,iv,"47")
        ig = ig_["23"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["36"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("47",ix_)
        arrset(pb.xnames,iv,"47")
        ig = ig_["56"]
        pbm.A[ig,iv] += Float64(-8.88)
        ig = ig_["57"]
        pbm.A[ig,iv] += Float64(-9.34)
        iv,ix_,_ = s2mpj_ii("47",ix_)
        arrset(pb.xnames,iv,"47")
        ig = ig_["58"]
        pbm.A[ig,iv] += Float64(-3.0)
        ig = ig_["59"]
        pbm.A[ig,iv] += Float64(12.0)
        iv,ix_,_ = s2mpj_ii("47",ix_)
        arrset(pb.xnames,iv,"47")
        ig = ig_["60"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["61"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("48",ix_)
        arrset(pb.xnames,iv,"48")
        ig = ig_["19"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["36"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("48",ix_)
        arrset(pb.xnames,iv,"48")
        ig = ig_["56"]
        pbm.A[ig,iv] += Float64(-9.47)
        ig = ig_["57"]
        pbm.A[ig,iv] += Float64(-9.46)
        iv,ix_,_ = s2mpj_ii("48",ix_)
        arrset(pb.xnames,iv,"48")
        ig = ig_["58"]
        pbm.A[ig,iv] += Float64(-3.0)
        ig = ig_["59"]
        pbm.A[ig,iv] += Float64(3.5)
        iv,ix_,_ = s2mpj_ii("48",ix_)
        arrset(pb.xnames,iv,"48")
        ig = ig_["60"]
        pbm.A[ig,iv] += Float64(.233)
        ig = ig_["61"]
        pbm.A[ig,iv] += Float64(-.358)
        iv,ix_,_ = s2mpj_ii("49",ix_)
        arrset(pb.xnames,iv,"49")
        ig = ig_["20"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["36"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("49",ix_)
        arrset(pb.xnames,iv,"49")
        ig = ig_["56"]
        pbm.A[ig,iv] += Float64(-9.07)
        ig = ig_["57"]
        pbm.A[ig,iv] += Float64(-9.21)
        iv,ix_,_ = s2mpj_ii("49",ix_)
        arrset(pb.xnames,iv,"49")
        ig = ig_["58"]
        pbm.A[ig,iv] += Float64(-3.0)
        ig = ig_["59"]
        pbm.A[ig,iv] += Float64(3.5)
        iv,ix_,_ = s2mpj_ii("49",ix_)
        arrset(pb.xnames,iv,"49")
        ig = ig_["60"]
        pbm.A[ig,iv] += Float64(.205)
        ig = ig_["61"]
        pbm.A[ig,iv] += Float64(-.333)
        iv,ix_,_ = s2mpj_ii("50",ix_)
        arrset(pb.xnames,iv,"50")
        ig = ig_["27"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["36"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("50",ix_)
        arrset(pb.xnames,iv,"50")
        ig = ig_["56"]
        pbm.A[ig,iv] += Float64(-9.27)
        ig = ig_["57"]
        pbm.A[ig,iv] += Float64(-9.14)
        iv,ix_,_ = s2mpj_ii("50",ix_)
        arrset(pb.xnames,iv,"50")
        ig = ig_["58"]
        pbm.A[ig,iv] += Float64(-3.0)
        ig = ig_["59"]
        pbm.A[ig,iv] += Float64(6.0)
        iv,ix_,_ = s2mpj_ii("50",ix_)
        arrset(pb.xnames,iv,"50")
        ig = ig_["60"]
        pbm.A[ig,iv] += Float64(.318)
        ig = ig_["61"]
        pbm.A[ig,iv] += Float64(-.509)
        iv,ix_,_ = s2mpj_ii("51",ix_)
        arrset(pb.xnames,iv,"51")
        ig = ig_["10"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["36"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("51",ix_)
        arrset(pb.xnames,iv,"51")
        ig = ig_["56"]
        pbm.A[ig,iv] += Float64(-9.78)
        ig = ig_["57"]
        pbm.A[ig,iv] += Float64(-9.79)
        iv,ix_,_ = s2mpj_ii("51",ix_)
        arrset(pb.xnames,iv,"51")
        ig = ig_["58"]
        pbm.A[ig,iv] += Float64(-3.0)
        ig = ig_["59"]
        pbm.A[ig,iv] += Float64(66.0)
        iv,ix_,_ = s2mpj_ii("51",ix_)
        arrset(pb.xnames,iv,"51")
        ig = ig_["60"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["61"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("52",ix_)
        arrset(pb.xnames,iv,"52")
        ig = ig_["56"]
        pbm.A[ig,iv] += Float64(-.426)
        ig = ig_["57"]
        pbm.A[ig,iv] += Float64(-.204)
        iv,ix_,_ = s2mpj_ii("52",ix_)
        arrset(pb.xnames,iv,"52")
        ig = ig_["58"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["C"]
        pbm.A[ig,iv] += Float64(.0924)
        iv,ix_,_ = s2mpj_ii("53",ix_)
        arrset(pb.xnames,iv,"53")
        ig = ig_["36"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["56"]
        pbm.A[ig,iv] += Float64(9.05)
        iv,ix_,_ = s2mpj_ii("53",ix_)
        arrset(pb.xnames,iv,"53")
        ig = ig_["57"]
        pbm.A[ig,iv] += Float64(9.05)
        ig = ig_["59"]
        pbm.A[ig,iv] += Float64(-9.5)
        iv,ix_,_ = s2mpj_ii("53",ix_)
        arrset(pb.xnames,iv,"53")
        ig = ig_["60"]
        pbm.A[ig,iv] += Float64(-.5)
        ig = ig_["61"]
        pbm.A[ig,iv] += Float64(.5)
        iv,ix_,_ = s2mpj_ii("53",ix_)
        arrset(pb.xnames,iv,"53")
        ig = ig_["73"]
        pbm.A[ig,iv] += Float64(-.36)
        ig = ig_["74"]
        pbm.A[ig,iv] += Float64(-.65)
        iv,ix_,_ = s2mpj_ii("53",ix_)
        arrset(pb.xnames,iv,"53")
        ig = ig_["C"]
        pbm.A[ig,iv] += Float64(-4.51)
        iv,ix_,_ = s2mpj_ii("54",ix_)
        arrset(pb.xnames,iv,"54")
        ig = ig_["9"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["26"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("55",ix_)
        arrset(pb.xnames,iv,"55")
        ig = ig_["9"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["37"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("56",ix_)
        arrset(pb.xnames,iv,"56")
        ig = ig_["10"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["37"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("57",ix_)
        arrset(pb.xnames,iv,"57")
        ig = ig_["37"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["C"]
        pbm.A[ig,iv] += Float64(-2.75)
        iv,ix_,_ = s2mpj_ii("58",ix_)
        arrset(pb.xnames,iv,"58")
        ig = ig_["11"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["38"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("58",ix_)
        arrset(pb.xnames,iv,"58")
        ig = ig_["63"]
        pbm.A[ig,iv] += Float64(-14.0)
        ig = ig_["64"]
        pbm.A[ig,iv] += Float64(14.0)
        iv,ix_,_ = s2mpj_ii("59",ix_)
        arrset(pb.xnames,iv,"59")
        ig = ig_["12"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["38"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("59",ix_)
        arrset(pb.xnames,iv,"59")
        ig = ig_["63"]
        pbm.A[ig,iv] += Float64(-.8)
        ig = ig_["64"]
        pbm.A[ig,iv] += Float64(.8)
        iv,ix_,_ = s2mpj_ii("60",ix_)
        arrset(pb.xnames,iv,"60")
        ig = ig_["38"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["63"]
        pbm.A[ig,iv] += Float64(2.0)
        iv,ix_,_ = s2mpj_ii("60",ix_)
        arrset(pb.xnames,iv,"60")
        ig = ig_["64"]
        pbm.A[ig,iv] += Float64(-3.0)
        ig = ig_["C"]
        pbm.A[ig,iv] += Float64(-4.2)
        iv,ix_,_ = s2mpj_ii("61",ix_)
        arrset(pb.xnames,iv,"61")
        ig = ig_["4"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["34"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("62",ix_)
        arrset(pb.xnames,iv,"62")
        ig = ig_["3"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["34"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("63",ix_)
        arrset(pb.xnames,iv,"63")
        ig = ig_["34"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["65"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("63",ix_)
        arrset(pb.xnames,iv,"63")
        ig = ig_["C"]
        pbm.A[ig,iv] += Float64(-3.6)
        iv,ix_,_ = s2mpj_ii("64",ix_)
        arrset(pb.xnames,iv,"64")
        ig = ig_["7"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["35"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("64",ix_)
        arrset(pb.xnames,iv,"64")
        ig = ig_["62"]
        pbm.A[ig,iv] += Float64(10.1)
        iv,ix_,_ = s2mpj_ii("65",ix_)
        arrset(pb.xnames,iv,"65")
        ig = ig_["8"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["35"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("65",ix_)
        arrset(pb.xnames,iv,"65")
        ig = ig_["62"]
        pbm.A[ig,iv] += Float64(12.63)
        iv,ix_,_ = s2mpj_ii("66",ix_)
        arrset(pb.xnames,iv,"66")
        ig = ig_["6"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["35"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("66",ix_)
        arrset(pb.xnames,iv,"66")
        ig = ig_["62"]
        pbm.A[ig,iv] += Float64(8.05)
        ig = ig_["66"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("67",ix_)
        arrset(pb.xnames,iv,"67")
        ig = ig_["5"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["35"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("67",ix_)
        arrset(pb.xnames,iv,"67")
        ig = ig_["62"]
        pbm.A[ig,iv] += Float64(6.9)
        ig = ig_["65"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("68",ix_)
        arrset(pb.xnames,iv,"68")
        ig = ig_["29"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["35"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("68",ix_)
        arrset(pb.xnames,iv,"68")
        ig = ig_["62"]
        pbm.A[ig,iv] += Float64(8.05)
        iv,ix_,_ = s2mpj_ii("69",ix_)
        arrset(pb.xnames,iv,"69")
        ig = ig_["28"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["35"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv,ix_,_ = s2mpj_ii("69",ix_)
        arrset(pb.xnames,iv,"69")
        ig = ig_["62"]
        pbm.A[ig,iv] += Float64(4.4)
        iv,ix_,_ = s2mpj_ii("70",ix_)
        arrset(pb.xnames,iv,"70")
        ig = ig_["35"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["62"]
        pbm.A[ig,iv] += Float64(-10.1)
        iv,ix_,_ = s2mpj_ii("70",ix_)
        arrset(pb.xnames,iv,"70")
        ig = ig_["C"]
        pbm.A[ig,iv] += Float64(-2.0)
        iv,ix_,_ = s2mpj_ii("71",ix_)
        arrset(pb.xnames,iv,"71")
        ig = ig_["39"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["41"]
        pbm.A[ig,iv] += Float64(-.325)
        iv,ix_,_ = s2mpj_ii("72",ix_)
        arrset(pb.xnames,iv,"72")
        ig = ig_["13"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["41"]
        pbm.A[ig,iv] += Float64(-4.153)
        iv,ix_,_ = s2mpj_ii("73",ix_)
        arrset(pb.xnames,iv,"73")
        ig = ig_["10"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["41"]
        pbm.A[ig,iv] += Float64(-4.316)
        iv,ix_,_ = s2mpj_ii("74",ix_)
        arrset(pb.xnames,iv,"74")
        ig = ig_["9"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["41"]
        pbm.A[ig,iv] += Float64(-3.814)
        iv,ix_,_ = s2mpj_ii("75",ix_)
        arrset(pb.xnames,iv,"75")
        ig = ig_["25"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["41"]
        pbm.A[ig,iv] += Float64(-3.808)
        iv,ix_,_ = s2mpj_ii("76",ix_)
        arrset(pb.xnames,iv,"76")
        ig = ig_["24"]
        pbm.A[ig,iv] += Float64(1.0)
        ig = ig_["41"]
        pbm.A[ig,iv] += Float64(-4.44)
        iv,ix_,_ = s2mpj_ii("77",ix_)
        arrset(pb.xnames,iv,"77")
        ig = ig_["40"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["41"]
        pbm.A[ig,iv] += Float64(1.42)
        iv,ix_,_ = s2mpj_ii("77",ix_)
        arrset(pb.xnames,iv,"77")
        ig = ig_["C"]
        pbm.A[ig,iv] += Float64(.04)
        iv,ix_,_ = s2mpj_ii("78",ix_)
        arrset(pb.xnames,iv,"78")
        ig = ig_["40"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("79",ix_)
        arrset(pb.xnames,iv,"79")
        ig = ig_["10"]
        pbm.A[ig,iv] += Float64(-.5)
        ig = ig_["13"]
        pbm.A[ig,iv] += Float64(-.5)
        iv,ix_,_ = s2mpj_ii("79",ix_)
        arrset(pb.xnames,iv,"79")
        ig = ig_["C"]
        pbm.A[ig,iv] += Float64(3.0)
        iv,ix_,_ = s2mpj_ii("80",ix_)
        arrset(pb.xnames,iv,"80")
        ig = ig_["41"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["C"]
        pbm.A[ig,iv] += Float64(.4)
        iv,ix_,_ = s2mpj_ii("81",ix_)
        arrset(pb.xnames,iv,"81")
        ig = ig_["41"]
        pbm.A[ig,iv] += Float64(1.0)
        iv,ix_,_ = s2mpj_ii("82",ix_)
        arrset(pb.xnames,iv,"82")
        ig = ig_["42"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["C"]
        pbm.A[ig,iv] += Float64(.0132)
        iv,ix_,_ = s2mpj_ii("83",ix_)
        arrset(pb.xnames,iv,"83")
        ig = ig_["43"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig = ig_["C"]
        pbm.A[ig,iv] += Float64(.01)
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
        pbm.gconst[ig_["65"]] = Float64(23.26)
        pbm.gconst[ig_["66"]] = Float64(5.25)
        pbm.gconst[ig_["67"]] = Float64(26.32)
        pbm.gconst[ig_["68"]] = Float64(21.05)
        pbm.gconst[ig_["69"]] = Float64(13.45)
        pbm.gconst[ig_["70"]] = Float64(2.58)
        pbm.gconst[ig_["71"]] = Float64(10.0)
        pbm.gconst[ig_["72"]] = Float64(10.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSQR", iet_)
        loaset(elftv,it,1,"X")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"D")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "E1"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-1.000000000))
        ename = "E2"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-0.865853667))
        ename = "E3"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-0.731707335))
        ename = "E4"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-0.597561002))
        ename = "E5"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-0.463414639))
        ename = "E6"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "6"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-0.329268306))
        ename = "E7"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "7"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-0.195121944))
        ename = "E8"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "8"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-0.060975611))
        ename = "E9"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "9"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.073170729))
        ename = "E10"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "10"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.207317069))
        ename = "E11"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "11"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.341463417))
        ename = "E12"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "12"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.475609750))
        ename = "E13"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "13"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.609756112))
        ename = "E14"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "14"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.743902445))
        ename = "E15"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "15"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.878048778))
        ename = "E16"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "16"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.012195110))
        ename = "E17"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "17"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.146341443))
        ename = "E18"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "18"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.280487776))
        ename = "E19"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "19"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.414634109))
        ename = "E20"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "20"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.548780441))
        ename = "E21"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "21"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.682926774))
        ename = "E22"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "22"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.817073226))
        ename = "E23"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "23"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.951219559))
        ename = "E24"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "24"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(2.085365772))
        ename = "E25"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "25"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(2.219512224))
        ename = "E26"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "26"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(2.353658438))
        ename = "E27"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "27"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(2.487804890))
        ename = "E28"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "28"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(2.621951103))
        ename = "E29"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "29"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(2.756097555))
        ename = "E30"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "30"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(2.890244007))
        ename = "E31"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "31"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(3.024390221))
        ename = "E32"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "32"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(3.158536673))
        ename = "E33"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "33"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(3.292682886))
        ename = "E34"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "34"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(3.426829338))
        ename = "E35"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "35"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(3.560975552))
        ename = "E36"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "36"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(3.695122004))
        ename = "E37"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "37"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(3.829268217))
        ename = "E38"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "38"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(3.963414669))
        ename = "E39"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "39"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(4.097560883))
        ename = "E40"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "40"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(4.231707096))
        ename = "E41"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "41"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(4.365853786))
        ename = "E42"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "42"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(4.500000000))
        ename = "E43"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "43"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(4.634146214))
        ename = "E44"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "44"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(4.768292904))
        ename = "E45"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "45"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(4.902439117))
        ename = "E46"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "46"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(5.036585331))
        ename = "E47"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "47"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(5.170731544))
        ename = "E48"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "48"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(5.304878235))
        ename = "E49"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "49"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(5.439024448))
        ename = "E50"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "50"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(5.573170662))
        ename = "E51"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "51"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(5.707316875))
        ename = "E52"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "52"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(5.841463566))
        ename = "E53"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "53"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(5.975609779))
        ename = "E54"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "54"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(6.109755993))
        ename = "E55"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "55"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(6.243902206))
        ename = "E56"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "56"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(6.378048897))
        ename = "E57"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "57"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(6.512195110))
        ename = "E58"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "58"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(6.646341324))
        ename = "E59"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "59"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(6.780488014))
        ename = "E60"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "60"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(6.914634228))
        ename = "E61"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "61"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(7.048780441))
        ename = "E62"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "62"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(7.182926655))
        ename = "E63"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "63"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(7.317073345))
        ename = "E64"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "64"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(7.451219559))
        ename = "E65"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "65"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(7.585365772))
        ename = "E66"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "66"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(7.719511986))
        ename = "E67"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "67"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(7.853658676))
        ename = "E68"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "68"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(7.987804890))
        ename = "E69"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "69"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(8.121951103))
        ename = "E70"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "70"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(8.256097794))
        ename = "E71"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "71"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(8.390243530))
        ename = "E72"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "72"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(8.524390221))
        ename = "E73"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "73"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(8.658536911))
        ename = "E74"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "74"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(8.792682648))
        ename = "E75"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "75"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(8.926829338))
        ename = "E76"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "76"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(9.060976028))
        ename = "E77"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "77"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(9.195121765))
        ename = "E78"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "78"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(9.329268456))
        ename = "E79"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "79"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(9.463414192))
        ename = "E80"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "80"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(9.597560883))
        ename = "E81"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "81"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(9.731707573))
        ename = "E82"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "82"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(9.865853310))
        ename = "E83"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "83"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="D",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(10.000000000))
        v_["1"] = 1
        v_["N"] = 83
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            ig = ig_["C"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
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
        pb.pbclass = "C-QLR2-MN-83-74"
        pb.x0          = zeros(Float64,pb.n)
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eSQR"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = pbm.elpar[iel_][1]*EV_[1]*EV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 2.0*pbm.elpar[iel_][1]*EV_[1]
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 2.0*pbm.elpar[iel_][1]
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

