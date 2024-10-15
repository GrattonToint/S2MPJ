function SSEBNLN(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : SSEBNLN
#    *********
#    The Power Generation problem for the SSGB.
# 
#    Source:
#    N. Gould, private communication.
# 
#    SIF input: Nick Gould, 23 October 1989
# 
#    classification = "C-LQR2-RN-194-96"
# 
#    period is the number of time periods
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "SSEBNLN"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["HOURS"] = 24
        v_["DAYS"] = 1
        v_["PERIOD"] = v_["HOURS"]*v_["DAYS"]
        v_["0"] = 0
        v_["1"] = 1
        v_["2"] = 2
        v_["Z1"] = 517.0
        v_["D1,1"] = 578.0
        v_["D1,2"] = 517.0
        v_["D1,3"] = 461.0
        v_["D1,4"] = 369.0
        v_["D1,5"] = 299.0
        v_["D1,6"] = 269.0
        v_["D1,7"] = 370.0
        v_["D1,8"] = 559.0
        v_["D1,9"] = 689.0
        v_["D1,10"] = 728.0
        v_["D1,11"] = 683.0
        v_["D1,12"] = 626.0
        v_["D1,13"] = 619.0
        v_["D1,14"] = 586.0
        v_["D1,15"] = 582.0
        v_["D1,16"] = 625.0
        v_["D1,17"] = 821.0
        v_["D1,18"] = 883.0
        v_["D1,19"] = 768.0
        v_["D1,20"] = 711.0
        v_["D1,21"] = 677.0
        v_["D1,22"] = 630.0
        v_["D1,23"] = 545.0
        v_["D1,24"] = 565.0
        v_["Z2"] = 400.0
        v_["D2,1"] = 631.0
        v_["D2,2"] = 574.0
        v_["D2,3"] = 521.0
        v_["D2,4"] = 446.0
        v_["D2,5"] = 359.0
        v_["D2,6"] = 336.0
        v_["D2,7"] = 420.0
        v_["D2,8"] = 588.0
        v_["D2,9"] = 697.0
        v_["D2,10"] = 732.0
        v_["D2,11"] = 713.0
        v_["D2,12"] = 682.0
        v_["D2,13"] = 695.0
        v_["D2,14"] = 651.0
        v_["D2,15"] = 645.0
        v_["D2,16"] = 664.0
        v_["D2,17"] = 816.0
        v_["D2,18"] = 858.0
        v_["D2,19"] = 760.0
        v_["D2,20"] = 700.0
        v_["D2,21"] = 659.0
        v_["D2,22"] = 623.0
        v_["D2,23"] = 517.0
        v_["D2,24"] = 542.0
        v_["Z3"] = 1017.0
        v_["D3,1"] = 582.0
        v_["D3,2"] = 501.0
        v_["D3,3"] = 443.0
        v_["D3,4"] = 367.0
        v_["D3,5"] = 288.0
        v_["D3,6"] = 265.0
        v_["D3,7"] = 349.0
        v_["D3,8"] = 503.0
        v_["D3,9"] = 663.0
        v_["D3,10"] = 651.0
        v_["D3,11"] = 625.0
        v_["D3,12"] = 596.0
        v_["D3,13"] = 608.0
        v_["D3,14"] = 566.0
        v_["D3,15"] = 555.0
        v_["D3,16"] = 584.0
        v_["D3,17"] = 763.0
        v_["D3,18"] = 803.0
        v_["D3,19"] = 710.0
        v_["D3,20"] = 648.0
        v_["D3,21"] = 626.0
        v_["D3,22"] = 590.0
        v_["D3,23"] = 486.0
        v_["D3,24"] = 540.0
        v_["Z4"] = 667.0
        v_["D4,1"] = 602.0
        v_["D4,2"] = 533.0
        v_["D4,3"] = 450.0
        v_["D4,4"] = 378.0
        v_["D4,5"] = 298.0
        v_["D4,6"] = 272.0
        v_["D4,7"] = 369.0
        v_["D4,8"] = 539.0
        v_["D4,9"] = 647.0
        v_["D4,10"] = 652.0
        v_["D4,11"] = 607.0
        v_["D4,12"] = 585.0
        v_["D4,13"] = 587.0
        v_["D4,14"] = 549.0
        v_["D4,15"] = 535.0
        v_["D4,16"] = 564.0
        v_["D4,17"] = 748.0
        v_["D4,18"] = 808.0
        v_["D4,19"] = 710.0
        v_["D4,20"] = 646.0
        v_["D4,21"] = 620.0
        v_["D4,22"] = 581.0
        v_["D4,23"] = 483.0
        v_["D4,24"] = 514.0
        v_["Z5"] = 600.0
        v_["D5,1"] = 579.0
        v_["D5,2"] = 518.0
        v_["D5,3"] = 447.0
        v_["D5,4"] = 355.0
        v_["D5,5"] = 284.0
        v_["D5,6"] = 261.0
        v_["D5,7"] = 348.0
        v_["D5,8"] = 530.0
        v_["D5,9"] = 644.0
        v_["D5,10"] = 648.0
        v_["D5,11"] = 607.0
        v_["D5,12"] = 570.0
        v_["D5,13"] = 577.0
        v_["D5,14"] = 536.0
        v_["D5,15"] = 544.0
        v_["D5,16"] = 554.0
        v_["D5,17"] = 716.0
        v_["D5,18"] = 765.0
        v_["D5,19"] = 676.0
        v_["D5,20"] = 631.0
        v_["D5,21"] = 576.0
        v_["D5,22"] = 528.0
        v_["D5,23"] = 445.0
        v_["D5,24"] = 520.0
        v_["Z6"] = 421.0
        v_["D6,1"] = 618.0
        v_["D6,2"] = 547.0
        v_["D6,3"] = 430.0
        v_["D6,4"] = 327.0
        v_["D6,5"] = 249.0
        v_["D6,6"] = 211.0
        v_["D6,7"] = 227.0
        v_["D6,8"] = 258.0
        v_["D6,9"] = 347.0
        v_["D6,10"] = 491.0
        v_["D6,11"] = 524.0
        v_["D6,12"] = 492.0
        v_["D6,13"] = 467.0
        v_["D6,14"] = 418.0
        v_["D6,15"] = 358.0
        v_["D6,16"] = 378.0
        v_["D6,17"] = 544.0
        v_["D6,18"] = 666.0
        v_["D6,19"] = 589.0
        v_["D6,20"] = 533.0
        v_["D6,21"] = 494.0
        v_["D6,22"] = 460.0
        v_["D6,23"] = 404.0
        v_["D6,24"] = 512.0
        v_["Z7"] = 425.0
        v_["D7,1"] = 615.0
        v_["D7,2"] = 587.0
        v_["D7,3"] = 450.0
        v_["D7,4"] = 320.0
        v_["D7,5"] = 235.0
        v_["D7,6"] = 198.0
        v_["D7,7"] = 195.0
        v_["D7,8"] = 173.0
        v_["D7,9"] = 197.0
        v_["D7,10"] = 349.0
        v_["D7,11"] = 441.0
        v_["D7,12"] = 459.0
        v_["D7,13"] = 485.0
        v_["D7,14"] = 445.0
        v_["D7,15"] = 410.0
        v_["D7,16"] = 421.0
        v_["D7,17"] = 568.0
        v_["D7,18"] = 643.0
        v_["D7,19"] = 596.0
        v_["D7,20"] = 566.0
        v_["D7,21"] = 541.0
        v_["D7,22"] = 532.0
        v_["D7,23"] = 454.0
        v_["D7,24"] = 511.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        iv,ix_,_  = (
              s2mpj_ii("V"*string(Int64(v_["0"]))*","*string(Int64(v_["HOURS"])),ix_))
        arrset(pb.xnames,iv,"V"*string(Int64(v_["0"]))*","*string(Int64(v_["HOURS"])))
        iv,ix_,_  = (
              s2mpj_ii("R"*string(Int64(v_["0"]))*","*string(Int64(v_["HOURS"])),ix_))
        arrset(pb.xnames,iv,"R"*string(Int64(v_["0"]))*","*string(Int64(v_["HOURS"])))
        for ID = Int64(v_["1"]):Int64(v_["DAYS"])
            for IH = Int64(v_["1"]):Int64(v_["HOURS"])
                iv,ix_,_ = s2mpj_ii("P1"*string(ID)*","*string(IH),ix_)
                arrset(pb.xnames,iv,"P1"*string(ID)*","*string(IH))
                iv,ix_,_ = s2mpj_ii("P2"*string(ID)*","*string(IH),ix_)
                arrset(pb.xnames,iv,"P2"*string(ID)*","*string(IH))
                iv,ix_,_ = s2mpj_ii("QH"*string(ID)*","*string(IH),ix_)
                arrset(pb.xnames,iv,"QH"*string(ID)*","*string(IH))
                iv,ix_,_ = s2mpj_ii("S"*string(ID)*","*string(IH),ix_)
                arrset(pb.xnames,iv,"S"*string(ID)*","*string(IH))
                iv,ix_,_ = s2mpj_ii("QG"*string(ID)*","*string(IH),ix_)
                arrset(pb.xnames,iv,"QG"*string(ID)*","*string(IH))
                iv,ix_,_ = s2mpj_ii("QP"*string(ID)*","*string(IH),ix_)
                arrset(pb.xnames,iv,"QP"*string(ID)*","*string(IH))
                iv,ix_,_ = s2mpj_ii("V"*string(ID)*","*string(IH),ix_)
                arrset(pb.xnames,iv,"V"*string(ID)*","*string(IH))
                iv,ix_,_ = s2mpj_ii("R"*string(ID)*","*string(IH),ix_)
                arrset(pb.xnames,iv,"R"*string(ID)*","*string(IH))
            end
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for ID = Int64(v_["1"]):Int64(v_["DAYS"])
            for IH = Int64(v_["1"]):Int64(v_["HOURS"])
                ig,ig_,_ = s2mpj_ii("OBJ",ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["P1"*string(ID)*","*string(IH)]
                pbm.A[ig,iv] += Float64(1000.0)
                iv = ix_["P2"*string(ID)*","*string(IH)]
                pbm.A[ig,iv] += Float64(1500.0)
                iv = ix_["QH"*string(ID)*","*string(IH)]
                pbm.A[ig,iv] += Float64(1200.0)
                iv = ix_["S"*string(ID)*","*string(IH)]
                pbm.A[ig,iv] += Float64(1200.0)
                iv = ix_["QG"*string(ID)*","*string(IH)]
                pbm.A[ig,iv] += Float64(1200.0)
                iv = ix_["QP"*string(ID)*","*string(IH)]
                pbm.A[ig,iv] += Float64(-1200.0)
            end
        end
        for ID = Int64(v_["1"]):Int64(v_["DAYS"])
            v_["P"] = -1+ID
            ig,ig_,_ = s2mpj_ii("H"*string(ID)*","*string(Int64(v_["1"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"H"*string(ID)*","*string(Int64(v_["1"])))
            iv = ix_["V"*string(ID)*","*string(Int64(v_["1"]))]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["V"*string(Int64(v_["P"]))*","*string(Int64(v_["HOURS"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("H"*string(ID)*","*string(Int64(v_["1"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"H"*string(ID)*","*string(Int64(v_["1"])))
            iv = ix_["S"*string(ID)*","*string(Int64(v_["1"]))]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["QH"*string(ID)*","*string(Int64(v_["1"]))]
            pbm.A[ig,iv] += Float64(1.0)
        end
        for ID = Int64(v_["1"]):Int64(v_["DAYS"])
            for IH = Int64(v_["2"]):Int64(v_["HOURS"])
                v_["IH-1"] = -1+IH
                ig,ig_,_ = s2mpj_ii("H"*string(ID)*","*string(IH),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"H"*string(ID)*","*string(IH))
                iv = ix_["V"*string(ID)*","*string(IH)]
                pbm.A[ig,iv] += Float64(1.0)
                iv = ix_["V"*string(ID)*","*string(Int64(v_["IH-1"]))]
                pbm.A[ig,iv] += Float64(-1.0)
                iv = ix_["S"*string(ID)*","*string(IH)]
                pbm.A[ig,iv] += Float64(1.0)
                iv = ix_["QH"*string(ID)*","*string(IH)]
                pbm.A[ig,iv] += Float64(1.0)
            end
        end
        for ID = Int64(v_["1"]):Int64(v_["DAYS"])
            v_["P"] = -1+ID
            ig,ig_,_ = s2mpj_ii("R"*string(ID)*","*string(Int64(v_["1"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"R"*string(ID)*","*string(Int64(v_["1"])))
            iv = ix_["R"*string(ID)*","*string(Int64(v_["1"]))]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["R"*string(Int64(v_["P"]))*","*string(Int64(v_["HOURS"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("R"*string(ID)*","*string(Int64(v_["1"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"R"*string(ID)*","*string(Int64(v_["1"])))
            iv = ix_["QG"*string(ID)*","*string(Int64(v_["1"]))]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["QP"*string(ID)*","*string(Int64(v_["1"]))]
            pbm.A[ig,iv] += Float64(-1.0)
        end
        for ID = Int64(v_["1"]):Int64(v_["DAYS"])
            for IH = Int64(v_["2"]):Int64(v_["HOURS"])
                v_["IH-1"] = -1+IH
                ig,ig_,_ = s2mpj_ii("R"*string(ID)*","*string(IH),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"R"*string(ID)*","*string(IH))
                iv = ix_["R"*string(ID)*","*string(IH)]
                pbm.A[ig,iv] += Float64(1.0)
                iv = ix_["R"*string(ID)*","*string(Int64(v_["IH-1"]))]
                pbm.A[ig,iv] += Float64(-1.0)
                iv = ix_["QG"*string(ID)*","*string(IH)]
                pbm.A[ig,iv] += Float64(1.0)
                iv = ix_["QP"*string(ID)*","*string(IH)]
                pbm.A[ig,iv] += Float64(-1.0)
            end
        end
        for ID = Int64(v_["1"]):Int64(v_["DAYS"])
            for IH = Int64(v_["1"]):Int64(v_["HOURS"])
                ig,ig_,_ = s2mpj_ii("D"*string(ID)*","*string(IH),ig_)
                arrset(gtype,ig,">=")
                arrset(pb.cnames,ig,"D"*string(ID)*","*string(IH))
                iv = ix_["P1"*string(ID)*","*string(IH)]
                pbm.A[ig,iv] += Float64(1.0)
                iv = ix_["P2"*string(ID)*","*string(IH)]
                pbm.A[ig,iv] += Float64(1.0)
                iv = ix_["QH"*string(ID)*","*string(IH)]
                pbm.A[ig,iv] += Float64(1.0)
                iv = ix_["QG"*string(ID)*","*string(IH)]
                pbm.A[ig,iv] += Float64(1.0)
                iv = ix_["QP"*string(ID)*","*string(IH)]
                pbm.A[ig,iv] += Float64(-1.33)
            end
        end
        for D = Int64(v_["1"]):Int64(v_["DAYS"])
            for H = Int64(v_["1"]):Int64(v_["HOURS"])
                ig,ig_,_ = s2mpj_ii("QG*QP"*string(D)*","*string(H),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"QG*QP"*string(D)*","*string(H))
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
        for ID = Int64(v_["1"]):Int64(v_["DAYS"])
            for IH = Int64(v_["1"]):Int64(v_["HOURS"])
                pbm.gconst[ig_["H"*string(ID)*","*string(IH)]] = Float64(v_["Z"*string(ID)])
            end
        end
        for ID = Int64(v_["1"]):Int64(v_["DAYS"])
            v_["0.01Z"] = 0.01*v_["Z"*string(ID)]
            for IH = Int64(v_["1"]):Int64(v_["HOURS"])
                pbm.gconst[ig_["R"*string(ID)*","*string(IH)]] = Float64(v_["0.01Z"])
            end
        end
        for ID = Int64(v_["1"]):Int64(v_["DAYS"])
            for IH = Int64(v_["1"]):Int64(v_["HOURS"])
                pbm.gconst[ig_["D"*string(ID)*","*string(IH)]]  = (
                      Float64(v_["D"*string(ID)*","*string(IH)]))
            end
        end
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower[ix_["V"*string(Int64(v_["0"]))*","*string(Int64(v_["HOURS"]))]]  = (
              240000.0)
        pb.xupper[ix_["V"*string(Int64(v_["0"]))*","*string(Int64(v_["HOURS"]))]]  = (
              240000.0)
        pb.xlower[ix_["R"*string(Int64(v_["0"]))*","*string(Int64(v_["HOURS"]))]]  = (
              3500.0)
        pb.xupper[ix_["R"*string(Int64(v_["0"]))*","*string(Int64(v_["HOURS"]))]]  = (
              3500.0)
        for ID = Int64(v_["1"]):Int64(v_["DAYS"])
            for IH = Int64(v_["1"]):Int64(v_["HOURS"])
                pb.xlower[ix_["P1"*string(ID)*","*string(IH)]] = 70.0
                pb.xlower[ix_["P2"*string(ID)*","*string(IH)]] = 90.0
                pb.xlower[ix_["QH"*string(ID)*","*string(IH)]] = 25.0
                pb.xlower[ix_["V"*string(ID)*","*string(IH)]] = 180000.0
                pb.xupper[ix_["P1"*string(ID)*","*string(IH)]] = 325.0
                pb.xupper[ix_["P2"*string(ID)*","*string(IH)]] = 290.0
                pb.xupper[ix_["QH"*string(ID)*","*string(IH)]] = 500.0
                pb.xupper[ix_["QP"*string(ID)*","*string(IH)]] = 225.0
                pb.xupper[ix_["QG"*string(ID)*","*string(IH)]] = 300.0
                pb.xupper[ix_["V"*string(ID)*","*string(IH)]] = 280000.0
                pb.xupper[ix_["R"*string(ID)*","*string(IH)]] = 6000.0
            end
        end
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        if haskey(ix_,"V"*string(Int64(v_["0"]))*","*string(Int64(v_["HOURS"])))
            pb.x0[ix_["V"*string(Int64(v_["0"]))*","*string(Int64(v_["HOURS"]))]]  = (
                  Float64(240000.0))
        else
            pb.y0[findfirst(x->x==ig_["V"*string(Int64(v_["0"]))*","*string(Int64(v_["HOURS"]))],pbm.congrps)] = Float64(240000.0)
        end
        if haskey(ix_,"R"*string(Int64(v_["0"]))*","*string(Int64(v_["HOURS"])))
            pb.x0[ix_["R"*string(Int64(v_["0"]))*","*string(Int64(v_["HOURS"]))]]  = (
                  Float64(3500.0))
        else
            pb.y0[findfirst(x->x==ig_["R"*string(Int64(v_["0"]))*","*string(Int64(v_["HOURS"]))],pbm.congrps)] = Float64(3500.0)
        end
        for ID = Int64(v_["1"]):Int64(v_["DAYS"])
            for IH = Int64(v_["1"]):Int64(v_["HOURS"])
                if haskey(ix_,"P1"*string(ID)*","*string(IH))
                    pb.x0[ix_["P1"*string(ID)*","*string(IH)]] = Float64(70.0)
                else
                    pb.y0[findfirst(x->x==ig_["P1"*string(ID)*","*string(IH)],pbm.congrps)]  = (
                          Float64(70.0))
                end
                if haskey(ix_,"P2"*string(ID)*","*string(IH))
                    pb.x0[ix_["P2"*string(ID)*","*string(IH)]] = Float64(90.0)
                else
                    pb.y0[findfirst(x->x==ig_["P2"*string(ID)*","*string(IH)],pbm.congrps)]  = (
                          Float64(90.0))
                end
                if haskey(ix_,"QH"*string(ID)*","*string(IH))
                    pb.x0[ix_["QH"*string(ID)*","*string(IH)]] = Float64(25.0)
                else
                    pb.y0[findfirst(x->x==ig_["QH"*string(ID)*","*string(IH)],pbm.congrps)]  = (
                          Float64(25.0))
                end
                if haskey(ix_,"QP"*string(ID)*","*string(IH))
                    pb.x0[ix_["QP"*string(ID)*","*string(IH)]] = Float64(225.0)
                else
                    pb.y0[findfirst(x->x==ig_["QP"*string(ID)*","*string(IH)],pbm.congrps)]  = (
                          Float64(225.0))
                end
                if haskey(ix_,"V"*string(ID)*","*string(IH))
                    pb.x0[ix_["V"*string(ID)*","*string(IH)]] = Float64(240000.0)
                else
                    pb.y0[findfirst(x->x==ig_["V"*string(ID)*","*string(IH)],pbm.congrps)]  = (
                          Float64(240000.0))
                end
                if haskey(ix_,"R"*string(ID)*","*string(IH))
                    pb.x0[ix_["R"*string(ID)*","*string(IH)]] = Float64(3500)
                else
                    pb.y0[findfirst(x->x==ig_["R"*string(ID)*","*string(IH)],pbm.congrps)]  = (
                          Float64(3500))
                end
            end
        end
        for ID = Int64(v_["1"]):Int64(v_["DAYS"])
            for IH = Int64(v_["1"]):Int64(v_["HOURS"])
            end
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "ePROD", iet_)
        loaset(elftv,it,1,"QP")
        loaset(elftv,it,2,"QG")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for ID = Int64(v_["1"]):Int64(v_["DAYS"])
            for IH = Int64(v_["1"]):Int64(v_["HOURS"])
                ename = "P"*string(ID)*","*string(IH)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"ePROD")
                arrset(ielftype,ie,iet_["ePROD"])
                vname = "QP"*string(ID)*","*string(IH)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="QP",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "QG"*string(ID)*","*string(IH)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="QG",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for D = Int64(v_["1"]):Int64(v_["DAYS"])
            for H = Int64(v_["1"]):Int64(v_["HOURS"])
                ig = ig_["QG*QP"*string(D)*","*string(H)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["P"*string(D)*","*string(H)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               1.617060D+07
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.clower[pb.nle+pb.neq+1:pb.m] = zeros(Float64,pb.nge)
        pb.cupper[1:pb.nge] = fill(Inf,pb.nge)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-LQR2-RN-194-96"
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

    elseif action == "ePROD"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[2]*EV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[2] = EV_[1]
            g_[1] = EV_[2]
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[2,1] = 1.0e+0
                H_[1,2] = H_[2,1]
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

