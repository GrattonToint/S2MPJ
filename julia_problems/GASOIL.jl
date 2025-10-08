function GASOIL(action::String,args::Union{Any}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : GASOIL
#    *********
# 
#    Determine the reaction coefficients for the catalytic cracking of gas oil
#    and other byproducts. The nonlinear model that describes the process is
# 
#      y_1' = - (theta_1 + theta_3 ) y_1^2
#      y_2' = theta_1 y_1^2 + theta_2 y_2
# 
#    with given initial conditions. The problem is to minimize
# 
#     sum{i=1,20} || y(tau_i,theta) - z_i||^2
# 
#    where the z_i are concentration measurements for y at times tau_i (i=1,20)
# 
#    This is problem 12 in the COPS (Version 2) collection of 
#    E. Dolan and J. More'
#    see "Benchmarking Optimization Software with COPS"
#    Argonne National Labs Technical Report ANL/MCS-246 (2000)
# 
#    SIF input: Nick Gould, November 2000
# 
#    classification = "C-COOR2-AN-V-V"
# 
#  The number of differential equations
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 8 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "GASOIL"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling GASOIL.")
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
        v_["NE"] = 2
        if nargin<1
            v_["NH"] = Int64(10);  #  SIF file default value
        else
            v_["NH"] = Int64(args[1]);
        end
#       Alternative values for the SIF file parameters:
# IE NH                  50             $-PARAMETER
# IE NH                  100            $-PARAMETER
# IE NH                  200            $-PARAMETER
# IE NH                  400            $-PARAMETER
        v_["NP"] = 3
        v_["NM"] = 21
        v_["NC"] = 4
        v_["RHO1"] = 0.0694318442
        v_["RHO2"] = 0.3300094782
        v_["RHO3"] = 0.6699905218
        v_["RHO4"] = 0.9305681558
        v_["TAU1"] = 0.0
        v_["TAU2"] = 0.025
        v_["TAU3"] = 0.05
        v_["TAU4"] = 0.075
        v_["TAU5"] = 0.10
        v_["TAU6"] = 0.125
        v_["TAU7"] = 0.150
        v_["TAU8"] = 0.175
        v_["TAU9"] = 0.20
        v_["TAU10"] = 0.225
        v_["TAU11"] = 0.250
        v_["TAU12"] = 0.30
        v_["TAU13"] = 0.35
        v_["TAU14"] = 0.40
        v_["TAU15"] = 0.45
        v_["TAU16"] = 0.50
        v_["TAU17"] = 0.55
        v_["TAU18"] = 0.65
        v_["TAU19"] = 0.75
        v_["TAU20"] = 0.85
        v_["TAU21"] = 0.95
        v_["TF"] = v_["TAU"*string(Int64(v_["NM"]))]
        v_["RNH"] = Float64(v_["NH"])
        v_["H"] = v_["TF"]/v_["RNH"]
        v_["Z1,1"] = 1.0000
        v_["Z1,2"] = 0.0000
        v_["Z2,1"] = 0.8105
        v_["Z2,2"] = 0.2000
        v_["Z3,1"] = 0.6208
        v_["Z3,2"] = 0.2886
        v_["Z4,1"] = 0.5258
        v_["Z4,2"] = 0.3010
        v_["Z5,1"] = 0.4345
        v_["Z5,2"] = 0.3215
        v_["Z6,1"] = 0.3903
        v_["Z6,2"] = 0.3123
        v_["Z7,1"] = 0.3342
        v_["Z7,2"] = 0.2716
        v_["Z8,1"] = 0.3034
        v_["Z8,2"] = 0.2551
        v_["Z9,1"] = 0.2735
        v_["Z9,2"] = 0.2258
        v_["Z10,1"] = 0.2405
        v_["Z10,2"] = 0.1959
        v_["Z11,1"] = 0.2283
        v_["Z11,2"] = 0.1789
        v_["Z12,1"] = 0.2071
        v_["Z12,2"] = 0.1457
        v_["Z13,1"] = 0.1669
        v_["Z13,2"] = 0.1198
        v_["Z14,1"] = 0.1530
        v_["Z14,2"] = 0.0909
        v_["Z15,1"] = 0.1339
        v_["Z15,2"] = 0.0719
        v_["Z16,1"] = 0.1265
        v_["Z16,2"] = 0.0561
        v_["Z17,1"] = 0.1200
        v_["Z17,2"] = 0.0460
        v_["Z18,1"] = 0.0990
        v_["Z18,2"] = 0.0280
        v_["Z19,1"] = 0.0870
        v_["Z19,2"] = 0.0190
        v_["Z20,1"] = 0.0770
        v_["Z20,2"] = 0.0140
        v_["Z21,1"] = 0.0690
        v_["Z21,2"] = 0.0100
        v_["BC1"] = 1.0
        v_["BC2"] = 0.0
        v_["1"] = 1
        v_["2"] = 2
        v_["3"] = 3
        v_["NH-1"] = -1+v_["NH"]
        v_["FACT0"] = 1.0
        for I = Int64(v_["1"]):Int64(v_["NC"])
            v_["RI"] = Float64(I)
            v_["I-1"] = -1+I
            v_["FACT"*string(I)] = v_["FACT"*string(Int64(v_["I-1"]))]*v_["RI"]
        end
        for I = Int64(v_["1"]):Int64(v_["NM"])
            v_["TAU/H"] = v_["TAU"*string(I)]/v_["H"]
            v_["IT/H"] = trunc(Int,v_["TAU/H"])
            v_["IT/H+1"] = 1+v_["IT/H"]
            v_["A"] = v_["IT/H+1"]
            v_["B"] = v_["NH"]
            v_["A"] = -1*v_["A"]
            v_["B"] = -1*v_["B"]
            v_["A"] = Float64(v_["A"])
            v_["ABSA"] = abs(v_["A"])
            v_["ABSA"] = trunc(Int,v_["ABSA"])
            v_["B"] = Float64(v_["B"])
            v_["ABSB"] = abs(v_["B"])
            v_["ABSB"] = trunc(Int,v_["ABSB"])
            v_["ABSA+B"] = v_["ABSA"]+v_["ABSB"]
            v_["A"] = v_["A"]+v_["ABSA+B"]
            v_["B"] = v_["B"]+v_["ABSA+B"]
            v_["A/B"] = trunc(Int,(v_["A"]/v_["B"]))
            v_["B/A"] = trunc(Int,(v_["B"]/v_["A"]))
            v_["SUM"] = v_["A/B"]+v_["B/A"]
            v_["A"] = v_["A"]*v_["A/B"]
            v_["B"] = v_["B"]*v_["B/A"]
            v_["MAXA,B"] = v_["A"]+v_["B"]
            v_["MAXA,B"] = trunc(Int,(v_["MAXA,B"]/v_["SUM"]))
            v_["MINA,B"] = v_["ABSA+B"]-v_["MAXA,B"]
            v_["ITAU"*string(I)] = Float64(v_["MINA,B"])
        end
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        irA   = Int64[]
        icA   = Int64[]
        valA  = Float64[]
        for I = Int64(v_["1"]):Int64(v_["NP"])
            iv,ix_,_ = s2mpj_ii("THETA"*string(I),ix_)
            arrset(pb.xnames,iv,"THETA"*string(I))
        end
        for I = Int64(v_["1"]):Int64(v_["NH"])
            for J = Int64(v_["1"]):Int64(v_["NE"])
                iv,ix_,_ = s2mpj_ii("V"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"V"*string(I)*","*string(J))
            end
            for K = Int64(v_["1"]):Int64(v_["NC"])
                for S = Int64(v_["1"]):Int64(v_["NE"])
                    iv,ix_,_ = s2mpj_ii("W"*string(I)*","*string(K)*","*string(S),ix_)
                    arrset(pb.xnames,iv,"W"*string(I)*","*string(K)*","*string(S))
                end
            end
            for J = Int64(v_["1"]):Int64(v_["NC"])
                for S = Int64(v_["1"]):Int64(v_["NE"])
                    iv,ix_,_ = s2mpj_ii("U"*string(I)*","*string(J)*","*string(S),ix_)
                    arrset(pb.xnames,iv,"U"*string(I)*","*string(J)*","*string(S))
                    iv,ix_,_ = s2mpj_ii("DU"*string(I)*","*string(J)*","*string(S),ix_)
                    arrset(pb.xnames,iv,"DU"*string(I)*","*string(J)*","*string(S))
                end
            end
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype = String[]
        for J = Int64(v_["1"]):Int64(v_["NM"])
            v_["RITAU"] = v_["ITAU"*string(J)]
            v_["I"] = trunc(Int,v_["RITAU"])
            v_["T"] = -1.0+v_["RITAU"]
            v_["T"] = v_["T"]*v_["H"]
            v_["DIFF"] = v_["TAU"*string(J)]-v_["T"]
            for S = Int64(v_["1"]):Int64(v_["NE"])
                v_["RATIO"] = v_["DIFF"]
                ig,ig_,_ = s2mpj_ii("OBJ"*string(J)*","*string(S),ig_)
                arrset(gtype,ig,"<>")
                push!(irA,ig)
                push!(icA,ix_["V"*string(Int64(v_["I"]))*","*string(S)])
                push!(valA,Float64(1.0))
                for K = Int64(v_["1"]):Int64(v_["NC"])
                    v_["COEF"] = v_["RATIO"]/v_["FACT"*string(K)]
                    ig,ig_,_ = s2mpj_ii("OBJ"*string(J)*","*string(S),ig_)
                    arrset(gtype,ig,"<>")
                    push!(irA,ig)
                    push!(icA,ix_["W"*string(Int64(v_["I"]))*","*string(K)*","*string(S)])
                    push!(valA,Float64(v_["COEF"]))
                    v_["RATIO"] = v_["RATIO"]*v_["DIFF"]
                    v_["RATIO"] = v_["RATIO"]/v_["H"]
                end
            end
        end
        for I = Int64(v_["1"]):Int64(v_["NH"])
            for J = Int64(v_["1"]):Int64(v_["NC"])
                v_["RH"] = v_["RHO"*string(J)]
                ig,ig_,_  = (
                      s2mpj_ii("U"*string(I)*","*string(J)*","*string(Int64(v_["1"])),ig_))
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"U"*string(I)*","*string(J)*","*string(Int64(v_["1"])))
                push!(irA,ig)
                push!(icA,ix_["U"*string(I)*","*string(J)*","*string(Int64(v_["1"]))])
                push!(valA,Float64(-1.0))
                push!(irA,ig)
                push!(icA,ix_["V"*string(I)*","*string(Int64(v_["1"]))])
                push!(valA,Float64(1.0))
                ig,ig_,_  = (
                      s2mpj_ii("U"*string(I)*","*string(J)*","*string(Int64(v_["2"])),ig_))
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"U"*string(I)*","*string(J)*","*string(Int64(v_["2"])))
                push!(irA,ig)
                push!(icA,ix_["U"*string(I)*","*string(J)*","*string(Int64(v_["2"]))])
                push!(valA,Float64(-1.0))
                push!(irA,ig)
                push!(icA,ix_["V"*string(I)*","*string(Int64(v_["2"]))])
                push!(valA,Float64(1.0))
                v_["PROD"] = v_["RH"]*v_["H"]
                for K = Int64(v_["1"]):Int64(v_["NC"])
                    v_["COEF"] = v_["PROD"]/v_["FACT"*string(K)]
                    ig,ig_,_  = (
                          s2mpj_ii("U"*string(I)*","*string(J)*","*string(Int64(v_["1"])),ig_))
                    arrset(gtype,ig,"==")
                    arrset(pb.cnames,ig,"U"*string(I)*","*string(J)*","*string(Int64(v_["1"])))
                    push!(irA,ig)
                    push!(icA,ix_["W"*string(I)*","*string(K)*","*string(Int64(v_["1"]))])
                    push!(valA,Float64(v_["COEF"]))
                    ig,ig_,_  = (
                          s2mpj_ii("U"*string(I)*","*string(J)*","*string(Int64(v_["2"])),ig_))
                    arrset(gtype,ig,"==")
                    arrset(pb.cnames,ig,"U"*string(I)*","*string(J)*","*string(Int64(v_["2"])))
                    push!(irA,ig)
                    push!(icA,ix_["W"*string(I)*","*string(K)*","*string(Int64(v_["2"]))])
                    push!(valA,Float64(v_["COEF"]))
                    v_["PROD"] = v_["PROD"]*v_["RH"]
                end
                ig,ig_,_  = (
                      s2mpj_ii("DU"*string(I)*","*string(J)*","*string(Int64(v_["1"])),ig_))
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"DU"*string(I)*","*string(J)*","*string(Int64(v_["1"])))
                push!(irA,ig)
                push!(icA,ix_["DU"*string(I)*","*string(J)*","*string(Int64(v_["1"]))])
                push!(valA,Float64(-1.0))
                ig,ig_,_  = (
                      s2mpj_ii("DU"*string(I)*","*string(J)*","*string(Int64(v_["2"])),ig_))
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"DU"*string(I)*","*string(J)*","*string(Int64(v_["2"])))
                push!(irA,ig)
                push!(icA,ix_["DU"*string(I)*","*string(J)*","*string(Int64(v_["2"]))])
                push!(valA,Float64(-1.0))
                v_["PROD"] = 1.0
                for K = Int64(v_["1"]):Int64(v_["NC"])
                    v_["K-1"] = -1+K
                    v_["COEF"] = v_["PROD"]/v_["FACT"*string(Int64(v_["K-1"]))]
                    ig,ig_,_  = (
                          s2mpj_ii("DU"*string(I)*","*string(J)*","*string(Int64(v_["1"])),ig_))
                    arrset(gtype,ig,"==")
                    arrset(pb.cnames,ig,"DU"*string(I)*","*string(J)*","*string(Int64(v_["1"])))
                    push!(irA,ig)
                    push!(icA,ix_["W"*string(I)*","*string(K)*","*string(Int64(v_["1"]))])
                    push!(valA,Float64(v_["COEF"]))
                    ig,ig_,_  = (
                          s2mpj_ii("DU"*string(I)*","*string(J)*","*string(Int64(v_["2"])),ig_))
                    arrset(gtype,ig,"==")
                    arrset(pb.cnames,ig,"DU"*string(I)*","*string(J)*","*string(Int64(v_["2"])))
                    push!(irA,ig)
                    push!(icA,ix_["W"*string(I)*","*string(K)*","*string(Int64(v_["2"]))])
                    push!(valA,Float64(v_["COEF"]))
                    v_["PROD"] = v_["PROD"]*v_["RH"]
                end
            end
        end
        for I = Int64(v_["1"]):Int64(v_["NH-1"])
            v_["I+1"] = 1+I
            for S = Int64(v_["1"]):Int64(v_["NE"])
                ig,ig_,_ = s2mpj_ii("C"*string(I)*","*string(S),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"C"*string(I)*","*string(S))
                push!(irA,ig)
                push!(icA,ix_["V"*string(I)*","*string(S)])
                push!(valA,Float64(1.0))
                push!(irA,ig)
                push!(icA,ix_["V"*string(Int64(v_["I+1"]))*","*string(S)])
                push!(valA,Float64(-1.0))
                for J = Int64(v_["1"]):Int64(v_["NC"])
                    v_["COEF"] = v_["H"]/v_["FACT"*string(J)]
                    ig,ig_,_ = s2mpj_ii("C"*string(I)*","*string(S),ig_)
                    arrset(gtype,ig,"==")
                    arrset(pb.cnames,ig,"C"*string(I)*","*string(S))
                    push!(irA,ig)
                    push!(icA,ix_["W"*string(I)*","*string(J)*","*string(S)])
                    push!(valA,Float64(v_["COEF"]))
                end
            end
        end
        for I = Int64(v_["1"]):Int64(v_["NH"])
            for J = Int64(v_["1"]):Int64(v_["NC"])
                for S = Int64(v_["1"]):Int64(v_["NE"])
                    ig,ig_,_ = s2mpj_ii("CO"*string(I)*","*string(J)*","*string(S),ig_)
                    arrset(gtype,ig,"==")
                    arrset(pb.cnames,ig,"CO"*string(I)*","*string(J)*","*string(S))
                    push!(irA,ig)
                    push!(icA,ix_["DU"*string(I)*","*string(J)*","*string(S)])
                    push!(valA,Float64(1.0))
                end
            end
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
        for J = Int64(v_["1"]):Int64(v_["NM"])
            for S = Int64(v_["1"]):Int64(v_["NE"])
                pbm.gconst[ig_["OBJ"*string(J)*","*string(S)]]  = (
                      Float64(v_["Z"*string(J)*","*string(S)]))
            end
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        for I = Int64(v_["1"]):Int64(v_["NP"])
            pb.xlower[ix_["THETA"*string(I)]] = 0.0
        end
        for S = Int64(v_["1"]):Int64(v_["NE"])
            pb.xlower[ix_["V"*string(Int64(v_["1"]))*","*string(S)]]  = (
                  v_["BC"*string(S)])
            pb.xupper[ix_["V"*string(Int64(v_["1"]))*","*string(S)]]  = (
                  v_["BC"*string(S)])
        end
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        for I = Int64(v_["1"]):Int64(v_["NP"])
            if haskey(ix_,"THETA"*string(I))
                pb.x0[ix_["THETA"*string(I)]] = Float64(0.0)
            else
                pb.y0[findfirst(x->x==ig_["THETA"*string(I)],pbm.congrps)] = Float64(0.0)
            end
        end
        for I = Int64(v_["1"]):Int64(v_["NH"])
            for J = Int64(v_["1"]):Int64(v_["NE"])
                if haskey(ix_,"V"*string(I)*","*string(J))
                    pb.x0[ix_["V"*string(I)*","*string(J)]] = Float64(0.0)
                else
                    pb.y0[findfirst(x->x==ig_["V"*string(I)*","*string(J)],pbm.congrps)]  = (
                          Float64(0.0))
                end
            end
        end
        v_["I1"] = 1
        v_["RITAU"] = v_["ITAU"*string(Int64(v_["1"]))]
        v_["I2"] = trunc(Int,v_["RITAU"])
        for I = Int64(v_["I1"]):Int64(v_["I2"])
            for S = Int64(v_["1"]):Int64(v_["NE"])
                if haskey(ix_,"V"*string(I)*","*string(S))
                    pb.x0[ix_["V"*string(I)*","*string(S)]] = Float64(v_["BC"*string(S)])
                else
                    pb.y0[findfirst(x->x==ig_["V"*string(I)*","*string(S)],pbm.congrps)]  = (
                          Float64(v_["BC"*string(S)]))
                end
                for J = Int64(v_["1"]):Int64(v_["NC"])
                    if haskey(ix_,"W"*string(I)*","*string(J)*","*string(S))
                        pb.x0[ix_["W"*string(I)*","*string(J)*","*string(S)]] = Float64(0.0)
                    else
                        pb.y0[findfirst(x->x==ig_["W"*string(I)*","*string(J)*","*string(S)],pbm.congrps)] = Float64(0.0)
                    end
                    if haskey(ix_,"U"*string(I)*","*string(J)*","*string(S))
                        pb.x0[ix_["U"*string(I)*","*string(J)*","*string(S)]]  = (
                              Float64(v_["BC"*string(S)]))
                    else
                        pb.y0[findfirst(x->x==ig_["U"*string(I)*","*string(J)*","*string(S)],pbm.congrps)] = Float64(v_["BC"*string(S)])
                    end
                    if haskey(ix_,"DU"*string(I)*","*string(J)*","*string(S))
                        pb.x0[ix_["DU"*string(I)*","*string(J)*","*string(S)]] = Float64(0.0)
                    else
                        pb.y0[findfirst(x->x==ig_["DU"*string(I)*","*string(J)*","*string(S)],pbm.congrps)] = Float64(0.0)
                    end
                end
            end
        end
        for K = Int64(v_["2"]):Int64(v_["NM"])
            v_["I1"] = 1+v_["I2"]
            v_["RITAU"] = v_["ITAU"*string(K)]
            v_["I2"] = trunc(Int,v_["RITAU"])
            for I = Int64(v_["I1"]):Int64(v_["I2"])
                v_["S"] = 1
                if haskey(ix_,"V"*string(I)*","*string(Int64(v_["S"])))
                    pb.x0[ix_["V"*string(I)*","*string(Int64(v_["S"]))]]  = (
                          Float64(v_["Z"*string(K)*","*string(Int64(v_["S"]))]))
                else
                    pb.y0[findfirst(x->x==ig_["V"*string(I)*","*string(Int64(v_["S"]))],pbm.congrps)] = Float64(v_["Z"*string(K)*","*string(Int64(v_["S"]))])
                end
                for J = Int64(v_["1"]):Int64(v_["NC"])
                    if haskey(ix_,"W"*string(I)*","*string(J)*","*string(Int64(v_["S"])))
                        pb.x0[ix_["W"*string(I)*","*string(J)*","*string(Int64(v_["S"]))]]  = (
                              Float64(0.0))
                    else
                        pb.y0[findfirst(x->x==ig_["W"*string(I)*","*string(J)*","*string(Int64(v_["S"]))],pbm.congrps)] = Float64(0.0)
                    end
                    if haskey(ix_,"U"*string(I)*","*string(J)*","*string(Int64(v_["S"])))
                        pb.x0[ix_["U"*string(I)*","*string(J)*","*string(Int64(v_["S"]))]]  = (
                              Float64(v_["Z"*string(K)*","*string(Int64(v_["S"]))]))
                    else
                        pb.y0[findfirst(x->x==ig_["U"*string(I)*","*string(J)*","*string(Int64(v_["S"]))],pbm.congrps)] = Float64(v_["Z"*string(K)*","*string(Int64(v_["S"]))])
                    end
                    if haskey(ix_,"DU"*string(I)*","*string(J)*","*string(Int64(v_["S"])))
                        pb.x0[ix_["DU"*string(I)*","*string(J)*","*string(Int64(v_["S"]))]]  = (
                              Float64(0.0))
                    else
                        pb.y0[findfirst(x->x==ig_["DU"*string(I)*","*string(J)*","*string(Int64(v_["S"]))],pbm.congrps)] = Float64(0.0)
                    end
                end
                v_["S"] = 2
                if haskey(ix_,"V"*string(I)*","*string(Int64(v_["S"])))
                    pb.x0[ix_["V"*string(I)*","*string(Int64(v_["S"]))]]  = (
                          Float64(v_["Z"*string(K)*","*string(Int64(v_["S"]))]))
                else
                    pb.y0[findfirst(x->x==ig_["V"*string(I)*","*string(Int64(v_["S"]))],pbm.congrps)] = Float64(v_["Z"*string(K)*","*string(Int64(v_["S"]))])
                end
                for J = Int64(v_["1"]):Int64(v_["NC"])
                    if haskey(ix_,"W"*string(I)*","*string(J)*","*string(Int64(v_["S"])))
                        pb.x0[ix_["W"*string(I)*","*string(J)*","*string(Int64(v_["S"]))]]  = (
                              Float64(0.0))
                    else
                        pb.y0[findfirst(x->x==ig_["W"*string(I)*","*string(J)*","*string(Int64(v_["S"]))],pbm.congrps)] = Float64(0.0)
                    end
                    if haskey(ix_,"U"*string(I)*","*string(J)*","*string(Int64(v_["S"])))
                        pb.x0[ix_["U"*string(I)*","*string(J)*","*string(Int64(v_["S"]))]]  = (
                              Float64(v_["Z"*string(K)*","*string(Int64(v_["S"]))]))
                    else
                        pb.y0[findfirst(x->x==ig_["U"*string(I)*","*string(J)*","*string(Int64(v_["S"]))],pbm.congrps)] = Float64(v_["Z"*string(K)*","*string(Int64(v_["S"]))])
                    end
                    if haskey(ix_,"DU"*string(I)*","*string(J)*","*string(Int64(v_["S"])))
                        pb.x0[ix_["DU"*string(I)*","*string(J)*","*string(Int64(v_["S"]))]]  = (
                              Float64(0.0))
                    else
                        pb.y0[findfirst(x->x==ig_["DU"*string(I)*","*string(J)*","*string(Int64(v_["S"]))],pbm.congrps)] = Float64(0.0)
                    end
                end
            end
        end
        v_["I1"] = 1+v_["I2"]
        v_["I2"] = v_["NH"]
        for I = Int64(v_["I1"]):Int64(v_["I2"])
            for S = Int64(v_["1"]):Int64(v_["NE"])
                if haskey(ix_,"V"*string(I)*","*string(S))
                    pb.x0[ix_["V"*string(I)*","*string(S)]]  = (
                          Float64(v_["Z"*string(Int64(v_["NM"]))*","*string(S)]))
                else
                    pb.y0[findfirst(x->x==ig_["V"*string(I)*","*string(S)],pbm.congrps)]  = (
                          Float64(v_["Z"*string(Int64(v_["NM"]))*","*string(S)]))
                end
                for J = Int64(v_["1"]):Int64(v_["NC"])
                    if haskey(ix_,"W"*string(I)*","*string(J)*","*string(S))
                        pb.x0[ix_["W"*string(I)*","*string(J)*","*string(S)]] = Float64(0.0)
                    else
                        pb.y0[findfirst(x->x==ig_["W"*string(I)*","*string(J)*","*string(S)],pbm.congrps)] = Float64(0.0)
                    end
                    if haskey(ix_,"U"*string(I)*","*string(J)*","*string(S))
                        pb.x0[ix_["U"*string(I)*","*string(J)*","*string(S)]]  = (
                              Float64(v_["Z"*string(Int64(v_["NM"]))*","*string(S)]))
                    else
                        pb.y0[findfirst(x->x==ig_["U"*string(I)*","*string(J)*","*string(S)],pbm.congrps)] = Float64(v_["Z"*string(Int64(v_["NM"]))*","*string(S)])
                    end
                    if haskey(ix_,"DU"*string(I)*","*string(J)*","*string(S))
                        pb.x0[ix_["DU"*string(I)*","*string(J)*","*string(S)]] = Float64(0.0)
                    else
                        pb.y0[findfirst(x->x==ig_["DU"*string(I)*","*string(J)*","*string(S)],pbm.congrps)] = Float64(0.0)
                    end
                end
            end
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "ePROD1", iet_)
        loaset(elftv,it,1,"THETA1")
        loaset(elftv,it,2,"THETA3")
        loaset(elftv,it,3,"U")
        it,iet_,_ = s2mpj_ii( "ePROD2", iet_)
        loaset(elftv,it,1,"THETA")
        loaset(elftv,it,2,"U")
        it,iet_,_ = s2mpj_ii( "ePROD3", iet_)
        loaset(elftv,it,1,"THETA")
        loaset(elftv,it,2,"U")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["NH"])
            for J = Int64(v_["1"]):Int64(v_["NC"])
                ename = "P1"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"ePROD1")
                arrset(ielftype,ie,iet_["ePROD1"])
                vname = "THETA"*string(Int64(v_["1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
                posev = findfirst(x->x=="THETA1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "THETA"*string(Int64(v_["3"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
                posev = findfirst(x->x=="THETA3",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "U"*string(I)*","*string(J)*","*string(Int64(v_["1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
                posev = findfirst(x->x=="U",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "P2"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"ePROD2")
                arrset(ielftype,ie,iet_["ePROD2"])
                vname = "THETA"*string(Int64(v_["1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
                posev = findfirst(x->x=="THETA",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "U"*string(I)*","*string(J)*","*string(Int64(v_["1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
                posev = findfirst(x->x=="U",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "P3"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"ePROD3")
                arrset(ielftype,ie,iet_["ePROD3"])
                vname = "THETA"*string(Int64(v_["2"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
                posev = findfirst(x->x=="THETA",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "U"*string(I)*","*string(J)*","*string(Int64(v_["2"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-Inf),Float64(Inf),nothing)
                posev = findfirst(x->x=="U",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gL2",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["NH"])
            for J = Int64(v_["1"]):Int64(v_["NC"])
                ig = ig_["CO"*string(I)*","*string(J)*","*string(Int64(v_["1"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["P1"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
                ig = ig_["CO"*string(I)*","*string(J)*","*string(Int64(v_["2"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["P2"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
                posel = posel+1
                loaset(pbm.grelt,ig,posel,ie_["P3"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel, 1.)
            end
        end
        for J = Int64(v_["1"]):Int64(v_["NM"])
            for S = Int64(v_["1"]):Int64(v_["NE"])
                ig = ig_["OBJ"*string(J)*","*string(S)]
                arrset(pbm.grftype,ig,"gL2")
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLUTION             5.23664D-03   $ (NH=50)
# LO SOLUTION             5.23659D-03   $ (NH=100)
# LO SOLUTION             5.23659D-03   $ (NH=200)
# LO SOLUTION             5.23659D-03   $ (NH=400)
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-COOR2-AN-V-V"
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

    elseif action == "ePROD1"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,2,3)
        IV_ =  zeros(Float64,2)
        U_[1,1] = U_[1,1]+1
        U_[1,2] = U_[1,2]+1
        U_[2,3] = U_[2,3]+1
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        f_   = IV_[1]*IV_[2]^2
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = IV_[2]^2
            g_[2] = 2.0*IV_[1]*IV_[2]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = 2.0*IV_[2]
                H_[2,1] = H_[1,2]
                H_[2,2] = 2.0*IV_[1]
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

    elseif action == "ePROD2"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = -EV_[1]*EV_[2]^2
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = -EV_[2]^2
            g_[2] = -2.0*EV_[1]*EV_[2]
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = -2.0*EV_[2]
                H_[2,1] = H_[1,2]
                H_[2,2] = -2.0*EV_[1]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "ePROD3"

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

    elseif action == "gL2"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_= GVAR_*GVAR_
        if nargout>1
            g_ = GVAR_+GVAR_
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = 2.0
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

