function TAX13322(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : TAX13322
#    --------
# 
#    An optimal income tax model with multidimensional taxpayer types,
#    due to Judd, Ma, Saunders & Su
# 
#    Source:
#    Kenneth L. Judd, Ma,  Michael A. Saunders and Che-Lin Su
#    "Optimal Income Taxation with Multidimensional Taxpayer Types"
#    Working Paper, Hoover Institute, Stanford University, 2017
# 
#    SIF input: Nick Gould, July 2018 based on the AMPL model pTAX5Dncl
# 
#    "If ever there was an example that exhibited the stupidity of SIF,
#     this is it. NIMG"
# 
#    classification = "C-OOR2-MN-72-1261"
# 
#    parameters
# 
#       Alternative values for the SIF file parameters:
# IE NA                  1              $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "TAX13322"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        if nargin<1
            v_["NA"] = Int64(1);  #  SIF file default value
        else
            v_["NA"] = Int64(args[1]);
        end
# IE NB                  3              $-PARAMETER
        if nargin<2
            v_["NB"] = Int64(3);  #  SIF file default value
        else
            v_["NB"] = Int64(args[2]);
        end
# IE NC                  3              $-PARAMETER
        if nargin<3
            v_["NC"] = Int64(3);  #  SIF file default value
        else
            v_["NC"] = Int64(args[3]);
        end
# IE ND                  2              $-PARAMETER
        if nargin<4
            v_["ND"] = Int64(2);  #  SIF file default value
        else
            v_["ND"] = Int64(args[4]);
        end
# IE NE                  2              $-PARAMETER
        if nargin<5
            v_["NE"] = Int64(2);  #  SIF file default value
        else
            v_["NE"] = Int64(args[5]);
        end
        v_["0"] = 0
        v_["1"] = 1
        v_["2"] = 2
        v_["3"] = 3
        v_["4"] = 4
        v_["5"] = 5
        v_["6"] = 6
        v_["ONE"] = 1.0e0
        v_["TWO"] = 2.0e0
        v_["THREE"] = 3.0e0
        v_["NBD"] = v_["NB"]*v_["ND"]
        v_["NCE"] = v_["NC"]*v_["NE"]
        v_["NP"] = v_["NBD"]*v_["NCE"]
        v_["NP"] = v_["NP"]*v_["NA"]
        v_["NPM1"] = -1+v_["NP"]
        v_["M"] = v_["NP"]*v_["NPM1"]
        v_["OMEGA1"] = v_["ONE"]/v_["TWO"]
        v_["OMEGA2"] = v_["TWO"]/v_["THREE"]
        v_["THETA1"] = v_["ONE"]/v_["THREE"]
        v_["THETA2"] = v_["ONE"]/v_["TWO"]
        v_["THETA3"] = v_["TWO"]/v_["THREE"]
        v_["PSI1"] = 1.0e0
        v_["PSI2"] = 1.5e0
        v_["W1"] = 2.0e0
        v_["W2"] = 2.5e0
        v_["W3"] = 3.0e0
        v_["W4"] = 3.5e0
        v_["W5"] = 4.0e0
        for I = Int64(v_["1"]):Int64(v_["NA"])
            for P = Int64(v_["1"]):Int64(v_["NBD"])
                for Q = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["LAM"*string(I)*","*string(P)*","*string(Q)] = 1.0e+0
                end
            end
        end
        v_["Q"] = v_["1"]
        v_["RA"*string(Int64(v_["Q"]))] = v_["ONE"]/v_["OMEGA1"]
        v_["Q"] = v_["2"]
        v_["RA"*string(Int64(v_["Q"]))] = v_["ONE"]/v_["OMEGA2"]
        v_["Q"] = v_["3"]
        v_["RA"*string(Int64(v_["Q"]))] = v_["ONE"]/v_["OMEGA1"]
        v_["Q"] = v_["4"]
        v_["RA"*string(Int64(v_["Q"]))] = v_["ONE"]/v_["OMEGA2"]
        v_["Q"] = v_["5"]
        v_["RA"*string(Int64(v_["Q"]))] = v_["ONE"]/v_["OMEGA1"]
        v_["Q"] = v_["6"]
        v_["RA"*string(Int64(v_["Q"]))] = v_["ONE"]/v_["OMEGA2"]
        for I = Int64(v_["1"]):Int64(v_["NA"])
            v_["LOGW"] = log(v_["W"*string(I)])
            v_["P"] = v_["1"]
            v_["-THETA"] = -1.0*v_["THETA1"]
            v_["RB"] = v_["LOGW"]*v_["-THETA"]
            v_["RE"] = exp(v_["RB"])
            v_["RB"] = v_["RE"]*v_["PSI1"]
            v_["RB"*string(I)*","*string(Int64(v_["P"]))] = v_["RB"]/v_["-THETA"]
            v_["P"] = v_["2"]
            v_["-THETA"] = -1.0*v_["THETA1"]
            v_["RB"] = v_["LOGW"]*v_["-THETA"]
            v_["RE"] = exp(v_["RB"])
            v_["RB"] = v_["RE"]*v_["PSI2"]
            v_["RB"*string(I)*","*string(Int64(v_["P"]))] = v_["RB"]/v_["-THETA"]
            v_["P"] = v_["3"]
            v_["-THETA"] = -1.0*v_["THETA2"]
            v_["RB"] = v_["LOGW"]*v_["-THETA"]
            v_["RE"] = exp(v_["RB"])
            v_["RB"] = v_["RE"]*v_["PSI1"]
            v_["RB"*string(I)*","*string(Int64(v_["P"]))] = v_["RB"]/v_["-THETA"]
            v_["P"] = v_["4"]
            v_["-THETA"] = -1.0*v_["THETA2"]
            v_["RB"] = v_["LOGW"]*v_["-THETA"]
            v_["RE"] = exp(v_["RB"])
            v_["RB"] = v_["RE"]*v_["PSI2"]
            v_["RB"*string(I)*","*string(Int64(v_["P"]))] = v_["RB"]/v_["-THETA"]
            v_["P"] = v_["5"]
            v_["-THETA"] = -1.0*v_["THETA3"]
            v_["RB"] = v_["LOGW"]*v_["-THETA"]
            v_["RE"] = exp(v_["RB"])
            v_["RB"] = v_["RE"]*v_["PSI1"]
            v_["RB"*string(I)*","*string(Int64(v_["P"]))] = v_["RB"]/v_["-THETA"]
            v_["P"] = v_["6"]
            v_["-THETA"] = -1.0*v_["THETA3"]
            v_["RB"] = v_["LOGW"]*v_["-THETA"]
            v_["RE"] = exp(v_["RB"])
            v_["RB"] = v_["RE"]*v_["PSI2"]
            v_["RB"*string(I)*","*string(Int64(v_["P"]))] = v_["RB"]/v_["-THETA"]
        end
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["NA"])
            for P = Int64(v_["1"]):Int64(v_["NBD"])
                for Q = Int64(v_["1"]):Int64(v_["NCE"])
                    iv,ix_,_ = s2mpj_ii("C"*string(I)*","*string(P)*","*string(Q),ix_)
                    arrset(pb.xnames,iv,"C"*string(I)*","*string(P)*","*string(Q))
                    iv,ix_,_ = s2mpj_ii("Y"*string(I)*","*string(P)*","*string(Q),ix_)
                    arrset(pb.xnames,iv,"Y"*string(I)*","*string(P)*","*string(Q))
                end
            end
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        for L = Int64(v_["1"]):Int64(v_["M"])
            ig,ig_,_ = s2mpj_ii("I"*string(L),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"I"*string(L))
        end
        for I = Int64(v_["1"]):Int64(v_["NA"])
            for P = Int64(v_["1"]):Int64(v_["NBD"])
                for Q = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["LAMBDA"] = v_["LAM"*string(I)*","*string(P)*","*string(Q)]
                    v_["-LAMBDA"] = -1.0e0*v_["LAMBDA"]
                    ig,ig_,_ = s2mpj_ii("T",ig_)
                    arrset(gtype,ig,">=")
                    arrset(pb.cnames,ig,"T")
                    iv = ix_["Y"*string(I)*","*string(P)*","*string(Q)]
                    pbm.A[ig,iv] += Float64(v_["LAMBDA"])
                    iv = ix_["C"*string(I)*","*string(P)*","*string(Q)]
                    pbm.A[ig,iv] += Float64(v_["-LAMBDA"])
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
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        for I = Int64(v_["1"]):Int64(v_["NA"])
            for P = Int64(v_["1"]):Int64(v_["NBD"])
                for Q = Int64(v_["1"]):Int64(v_["NCE"])
                    pb.xlower[ix_["C"*string(I)*","*string(P)*","*string(Q)]] = 0.1e0
                    pb.xlower[ix_["Y"*string(I)*","*string(P)*","*string(Q)]] = 0.1e0
                end
            end
        end
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(0.1e0),pb.n)
        pb.y0 = fill(Float64(0.1e0),pb.m)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eA1", iet_)
        loaset(elftv,it,1,"C")
        it,iet_,_ = s2mpj_ii( "eA2", iet_)
        loaset(elftv,it,1,"C")
        it,iet_,_ = s2mpj_ii( "eA3", iet_)
        loaset(elftv,it,1,"C")
        it,iet_,_ = s2mpj_ii( "eA4", iet_)
        loaset(elftv,it,1,"C")
        it,iet_,_ = s2mpj_ii( "eA5", iet_)
        loaset(elftv,it,1,"C")
        it,iet_,_ = s2mpj_ii( "eA6", iet_)
        loaset(elftv,it,1,"C")
        it,iet_,_ = s2mpj_ii( "eB1", iet_)
        loaset(elftv,it,1,"Y")
        it,iet_,_ = s2mpj_ii( "eB2", iet_)
        loaset(elftv,it,1,"Y")
        it,iet_,_ = s2mpj_ii( "eB3", iet_)
        loaset(elftv,it,1,"Y")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["NA"])
            for P = Int64(v_["1"]):Int64(v_["NBD"])
                for Q = Int64(v_["1"]):Int64(v_["NCE"])
                    ename = "A1-"*string(I)*","*string(P)*","*string(Q)
                    ie,ie_,_  = s2mpj_ii(ename,ie_)
                    arrset(pbm.elftype,ie,"eA1")
                    arrset(ielftype,ie,iet_["eA1"])
                    vname = "C"*string(I)*","*string(P)*","*string(Q)
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.1e0))
                    posev = findfirst(x->x=="C",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                end
            end
        end
        for I = Int64(v_["1"]):Int64(v_["NA"])
            for P = Int64(v_["1"]):Int64(v_["NBD"])
                for Q = Int64(v_["1"]):Int64(v_["NCE"])
                    ename = "A2-"*string(I)*","*string(P)*","*string(Q)
                    ie,ie_,_  = s2mpj_ii(ename,ie_)
                    arrset(pbm.elftype,ie,"eA2")
                    arrset(ielftype,ie,iet_["eA2"])
                    vname = "C"*string(I)*","*string(P)*","*string(Q)
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.1e0))
                    posev = findfirst(x->x=="C",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                end
            end
        end
        for I = Int64(v_["1"]):Int64(v_["NA"])
            for P = Int64(v_["1"]):Int64(v_["NBD"])
                for Q = Int64(v_["1"]):Int64(v_["NCE"])
                    ename = "A3-"*string(I)*","*string(P)*","*string(Q)
                    ie,ie_,_  = s2mpj_ii(ename,ie_)
                    arrset(pbm.elftype,ie,"eA3")
                    arrset(ielftype,ie,iet_["eA3"])
                    vname = "C"*string(I)*","*string(P)*","*string(Q)
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.1e0))
                    posev = findfirst(x->x=="C",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                end
            end
        end
        for I = Int64(v_["1"]):Int64(v_["NA"])
            for P = Int64(v_["1"]):Int64(v_["NBD"])
                for Q = Int64(v_["1"]):Int64(v_["NCE"])
                    ename = "A4-"*string(I)*","*string(P)*","*string(Q)
                    ie,ie_,_  = s2mpj_ii(ename,ie_)
                    arrset(pbm.elftype,ie,"eA4")
                    arrset(ielftype,ie,iet_["eA4"])
                    vname = "C"*string(I)*","*string(P)*","*string(Q)
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.1e0))
                    posev = findfirst(x->x=="C",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                end
            end
        end
        for I = Int64(v_["1"]):Int64(v_["NA"])
            for P = Int64(v_["1"]):Int64(v_["NBD"])
                for Q = Int64(v_["1"]):Int64(v_["NCE"])
                    ename = "A5-"*string(I)*","*string(P)*","*string(Q)
                    ie,ie_,_  = s2mpj_ii(ename,ie_)
                    arrset(pbm.elftype,ie,"eA5")
                    arrset(ielftype,ie,iet_["eA5"])
                    vname = "C"*string(I)*","*string(P)*","*string(Q)
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.1e0))
                    posev = findfirst(x->x=="C",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                end
            end
        end
        for I = Int64(v_["1"]):Int64(v_["NA"])
            for P = Int64(v_["1"]):Int64(v_["NBD"])
                for Q = Int64(v_["1"]):Int64(v_["NCE"])
                    ename = "A6-"*string(I)*","*string(P)*","*string(Q)
                    ie,ie_,_  = s2mpj_ii(ename,ie_)
                    arrset(pbm.elftype,ie,"eA6")
                    arrset(ielftype,ie,iet_["eA6"])
                    vname = "C"*string(I)*","*string(P)*","*string(Q)
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.1e0))
                    posev = findfirst(x->x=="C",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                end
            end
        end
        for I = Int64(v_["1"]):Int64(v_["NA"])
            for P = Int64(v_["1"]):Int64(v_["NBD"])
                for Q = Int64(v_["1"]):Int64(v_["NCE"])
                    ename = "B1-"*string(I)*","*string(P)*","*string(Q)
                    ie,ie_,_  = s2mpj_ii(ename,ie_)
                    arrset(pbm.elftype,ie,"eB1")
                    arrset(ielftype,ie,iet_["eB1"])
                    vname = "Y"*string(I)*","*string(P)*","*string(Q)
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.1e0))
                    posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                end
            end
        end
        for I = Int64(v_["1"]):Int64(v_["NA"])
            for P = Int64(v_["1"]):Int64(v_["NBD"])
                for Q = Int64(v_["1"]):Int64(v_["NCE"])
                    ename = "B2-"*string(I)*","*string(P)*","*string(Q)
                    ie,ie_,_  = s2mpj_ii(ename,ie_)
                    arrset(pbm.elftype,ie,"eB2")
                    arrset(ielftype,ie,iet_["eB2"])
                    vname = "Y"*string(I)*","*string(P)*","*string(Q)
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.1e0))
                    posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                end
            end
        end
        for I = Int64(v_["1"]):Int64(v_["NA"])
            for P = Int64(v_["1"]):Int64(v_["NBD"])
                for Q = Int64(v_["1"]):Int64(v_["NCE"])
                    ename = "B3-"*string(I)*","*string(P)*","*string(Q)
                    ie,ie_,_  = s2mpj_ii(ename,ie_)
                    arrset(pbm.elftype,ie,"eB3")
                    arrset(ielftype,ie,iet_["eB3"])
                    vname = "Y"*string(I)*","*string(P)*","*string(Q)
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.1e0))
                    posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                end
            end
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["NA"])
            for P = Int64(v_["1"]):Int64(v_["NBD"])
                v_["Q"] = v_["1"]
                v_["R"]  = (
                      v_["RA"*string(Int64(v_["Q"]))]*v_["LAM"*string(I)*","*string(P)*","*string(Int64(v_["Q"]))])
                ig = ig_["OBJ"]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(P)*","*string(Int64(v_["Q"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["R"]))
                v_["Q"] = v_["2"]
                v_["R"]  = (
                      v_["RA"*string(Int64(v_["Q"]))]*v_["LAM"*string(I)*","*string(P)*","*string(Int64(v_["Q"]))])
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(P)*","*string(Int64(v_["Q"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["R"]))
                v_["Q"] = v_["3"]
                v_["R"]  = (
                      v_["RA"*string(Int64(v_["Q"]))]*v_["LAM"*string(I)*","*string(P)*","*string(Int64(v_["Q"]))])
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(P)*","*string(Int64(v_["Q"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["R"]))
                v_["Q"] = v_["4"]
                v_["R"]  = (
                      v_["RA"*string(Int64(v_["Q"]))]*v_["LAM"*string(I)*","*string(P)*","*string(Int64(v_["Q"]))])
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(P)*","*string(Int64(v_["Q"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["R"]))
                v_["Q"] = v_["5"]
                v_["R"]  = (
                      v_["RA"*string(Int64(v_["Q"]))]*v_["LAM"*string(I)*","*string(P)*","*string(Int64(v_["Q"]))])
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(P)*","*string(Int64(v_["Q"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["R"]))
                v_["Q"] = v_["6"]
                v_["R"]  = (
                      v_["RA"*string(Int64(v_["Q"]))]*v_["LAM"*string(I)*","*string(P)*","*string(Int64(v_["Q"]))])
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(P)*","*string(Int64(v_["Q"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["R"]))
            end
        end
        for I = Int64(v_["1"]):Int64(v_["NA"])
            v_["LOGW"] = log(v_["W"*string(I)])
            for Q = Int64(v_["1"]):Int64(v_["NCE"])
                v_["P"] = v_["1"]
                v_["R"]  = (
                      v_["RB"*string(I)*","*string(Int64(v_["P"]))]*v_["LAM"*string(I)*","*string(Int64(v_["P"]))*","*string(Q)])
                ig = ig_["OBJ"]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Q)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["R"]))
                v_["P"] = v_["2"]
                v_["R"]  = (
                      v_["RB"*string(I)*","*string(Int64(v_["P"]))]*v_["LAM"*string(I)*","*string(Int64(v_["P"]))*","*string(Q)])
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Q)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["R"]))
                v_["P"] = v_["3"]
                v_["R"]  = (
                      v_["RB"*string(I)*","*string(Int64(v_["P"]))]*v_["LAM"*string(I)*","*string(Int64(v_["P"]))*","*string(Q)])
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Q)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["R"]))
                v_["P"] = v_["4"]
                v_["R"]  = (
                      v_["RB"*string(I)*","*string(Int64(v_["P"]))]*v_["LAM"*string(I)*","*string(Int64(v_["P"]))*","*string(Q)])
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Q)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["R"]))
                v_["P"] = v_["5"]
                v_["R"]  = (
                      v_["RB"*string(I)*","*string(Int64(v_["P"]))]*v_["LAM"*string(I)*","*string(Int64(v_["P"]))*","*string(Q)])
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Q)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["R"]))
                v_["P"] = v_["6"]
                v_["R"]  = (
                      v_["RB"*string(I)*","*string(Int64(v_["P"]))]*v_["LAM"*string(I)*","*string(Int64(v_["P"]))*","*string(Q)])
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Q)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["R"]))
            end
        end
        v_["L"] = v_["0"]
        for I = Int64(v_["1"]):Int64(v_["NA"])
            v_["I+1"] = 1+I
            v_["I-1"] = -1+I
            v_["LOGW"] = log(v_["W"*string(I)])
            v_["P"] = v_["1"]
            v_["P+1"] = 1+v_["P"]
            v_["P-1"] = -1+v_["P"]
            v_["RB"] = v_["RB"*string(I)*","*string(Int64(v_["P"]))]
            v_["-RB"] = -1.0*v_["RB"]
            v_["RA"] = v_["RA"*string(Int64(v_["1"]))]
            v_["-RA"] = -1.0*v_["RA"]
            for S = Int64(v_["1"]):Int64(v_["P-1"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for S = Int64(v_["P+1"]):Int64(v_["NBD"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["1"]):Int64(v_["I-1"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["I+1"]):Int64(v_["NA"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for T = Int64(v_["2"]):Int64(v_["NCE"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            v_["RA"] = v_["RA"*string(Int64(v_["2"]))]
            v_["-RA"] = -1.0*v_["RA"]
            for S = Int64(v_["1"]):Int64(v_["P-1"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for S = Int64(v_["P+1"]):Int64(v_["NBD"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["1"]):Int64(v_["I-1"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["I+1"]):Int64(v_["NA"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for T = Int64(v_["1"]):Int64(v_["1"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            for T = Int64(v_["3"]):Int64(v_["NCE"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            v_["RA"] = v_["RA"*string(Int64(v_["3"]))]
            v_["-RA"] = -1.0*v_["RA"]
            for S = Int64(v_["1"]):Int64(v_["P-1"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for S = Int64(v_["P+1"]):Int64(v_["NBD"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["1"]):Int64(v_["I-1"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["I+1"]):Int64(v_["NA"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for T = Int64(v_["1"]):Int64(v_["2"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            for T = Int64(v_["4"]):Int64(v_["NCE"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            v_["RA"] = v_["RA"*string(Int64(v_["4"]))]
            v_["-RA"] = -1.0*v_["RA"]
            for S = Int64(v_["1"]):Int64(v_["P-1"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for S = Int64(v_["P+1"]):Int64(v_["NBD"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["1"]):Int64(v_["I-1"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["I+1"]):Int64(v_["NA"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for T = Int64(v_["1"]):Int64(v_["3"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            for T = Int64(v_["5"]):Int64(v_["NCE"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            v_["RA"] = v_["RA"*string(Int64(v_["5"]))]
            v_["-RA"] = -1.0*v_["RA"]
            for S = Int64(v_["1"]):Int64(v_["P-1"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for S = Int64(v_["P+1"]):Int64(v_["NBD"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["1"]):Int64(v_["I-1"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["I+1"]):Int64(v_["NA"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for T = Int64(v_["1"]):Int64(v_["4"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            for T = Int64(v_["6"]):Int64(v_["NCE"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            v_["RA"] = v_["RA"*string(Int64(v_["6"]))]
            v_["-RA"] = -1.0*v_["RA"]
            for S = Int64(v_["1"]):Int64(v_["P-1"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for S = Int64(v_["P+1"]):Int64(v_["NBD"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["1"]):Int64(v_["I-1"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["I+1"]):Int64(v_["NA"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for T = Int64(v_["1"]):Int64(v_["5"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            v_["P"] = v_["2"]
            v_["P+1"] = 1+v_["P"]
            v_["P-1"] = -1+v_["P"]
            v_["RB"] = v_["RB"*string(I)*","*string(Int64(v_["P"]))]
            v_["-RB"] = -1.0*v_["RB"]
            v_["RA"] = v_["RA"*string(Int64(v_["1"]))]
            v_["-RA"] = -1.0*v_["RA"]
            for S = Int64(v_["1"]):Int64(v_["P-1"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for S = Int64(v_["P+1"]):Int64(v_["NBD"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["1"]):Int64(v_["I-1"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["I+1"]):Int64(v_["NA"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for T = Int64(v_["2"]):Int64(v_["NCE"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            v_["RA"] = v_["RA"*string(Int64(v_["2"]))]
            v_["-RA"] = -1.0*v_["RA"]
            for S = Int64(v_["1"]):Int64(v_["P-1"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for S = Int64(v_["P+1"]):Int64(v_["NBD"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["1"]):Int64(v_["I-1"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["I+1"]):Int64(v_["NA"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for T = Int64(v_["1"]):Int64(v_["1"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            for T = Int64(v_["3"]):Int64(v_["NCE"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            v_["RA"] = v_["RA"*string(Int64(v_["3"]))]
            v_["-RA"] = -1.0*v_["RA"]
            for S = Int64(v_["1"]):Int64(v_["P-1"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for S = Int64(v_["P+1"]):Int64(v_["NBD"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["1"]):Int64(v_["I-1"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["I+1"]):Int64(v_["NA"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for T = Int64(v_["1"]):Int64(v_["2"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            for T = Int64(v_["4"]):Int64(v_["NCE"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            v_["RA"] = v_["RA"*string(Int64(v_["4"]))]
            v_["-RA"] = -1.0*v_["RA"]
            for S = Int64(v_["1"]):Int64(v_["P-1"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for S = Int64(v_["P+1"]):Int64(v_["NBD"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["1"]):Int64(v_["I-1"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["I+1"]):Int64(v_["NA"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for T = Int64(v_["1"]):Int64(v_["3"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            for T = Int64(v_["5"]):Int64(v_["NCE"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            v_["RA"] = v_["RA"*string(Int64(v_["5"]))]
            v_["-RA"] = -1.0*v_["RA"]
            for S = Int64(v_["1"]):Int64(v_["P-1"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for S = Int64(v_["P+1"]):Int64(v_["NBD"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["1"]):Int64(v_["I-1"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["I+1"]):Int64(v_["NA"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for T = Int64(v_["1"]):Int64(v_["4"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            for T = Int64(v_["6"]):Int64(v_["NCE"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            v_["RA"] = v_["RA"*string(Int64(v_["6"]))]
            v_["-RA"] = -1.0*v_["RA"]
            for S = Int64(v_["1"]):Int64(v_["P-1"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for S = Int64(v_["P+1"]):Int64(v_["NBD"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["1"]):Int64(v_["I-1"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["I+1"]):Int64(v_["NA"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B1-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for T = Int64(v_["1"]):Int64(v_["5"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B1-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            v_["P"] = v_["3"]
            v_["P+1"] = 1+v_["P"]
            v_["P-1"] = -1+v_["P"]
            v_["RB"] = v_["RB"*string(I)*","*string(Int64(v_["P"]))]
            v_["-RB"] = -1.0*v_["RB"]
            v_["RA"] = v_["RA"*string(Int64(v_["1"]))]
            v_["-RA"] = -1.0*v_["RA"]
            for S = Int64(v_["1"]):Int64(v_["P-1"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for S = Int64(v_["P+1"]):Int64(v_["NBD"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["1"]):Int64(v_["I-1"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["I+1"]):Int64(v_["NA"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for T = Int64(v_["2"]):Int64(v_["NCE"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            v_["RA"] = v_["RA"*string(Int64(v_["2"]))]
            v_["-RA"] = -1.0*v_["RA"]
            for S = Int64(v_["1"]):Int64(v_["P-1"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for S = Int64(v_["P+1"]):Int64(v_["NBD"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["1"]):Int64(v_["I-1"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["I+1"]):Int64(v_["NA"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for T = Int64(v_["1"]):Int64(v_["1"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            for T = Int64(v_["3"]):Int64(v_["NCE"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            v_["RA"] = v_["RA"*string(Int64(v_["3"]))]
            v_["-RA"] = -1.0*v_["RA"]
            for S = Int64(v_["1"]):Int64(v_["P-1"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for S = Int64(v_["P+1"]):Int64(v_["NBD"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["1"]):Int64(v_["I-1"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["I+1"]):Int64(v_["NA"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for T = Int64(v_["1"]):Int64(v_["2"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            for T = Int64(v_["4"]):Int64(v_["NCE"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            v_["RA"] = v_["RA"*string(Int64(v_["4"]))]
            v_["-RA"] = -1.0*v_["RA"]
            for S = Int64(v_["1"]):Int64(v_["P-1"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for S = Int64(v_["P+1"]):Int64(v_["NBD"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["1"]):Int64(v_["I-1"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["I+1"]):Int64(v_["NA"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for T = Int64(v_["1"]):Int64(v_["3"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            for T = Int64(v_["5"]):Int64(v_["NCE"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            v_["RA"] = v_["RA"*string(Int64(v_["5"]))]
            v_["-RA"] = -1.0*v_["RA"]
            for S = Int64(v_["1"]):Int64(v_["P-1"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for S = Int64(v_["P+1"]):Int64(v_["NBD"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["1"]):Int64(v_["I-1"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["I+1"]):Int64(v_["NA"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for T = Int64(v_["1"]):Int64(v_["4"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            for T = Int64(v_["6"]):Int64(v_["NCE"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            v_["RA"] = v_["RA"*string(Int64(v_["6"]))]
            v_["-RA"] = -1.0*v_["RA"]
            for S = Int64(v_["1"]):Int64(v_["P-1"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for S = Int64(v_["P+1"]):Int64(v_["NBD"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["1"]):Int64(v_["I-1"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["I+1"]):Int64(v_["NA"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for T = Int64(v_["1"]):Int64(v_["5"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            v_["P"] = v_["4"]
            v_["P+1"] = 1+v_["P"]
            v_["P-1"] = -1+v_["P"]
            v_["RB"] = v_["RB"*string(I)*","*string(Int64(v_["P"]))]
            v_["-RB"] = -1.0*v_["RB"]
            v_["RA"] = v_["RA"*string(Int64(v_["1"]))]
            v_["-RA"] = -1.0*v_["RA"]
            for S = Int64(v_["1"]):Int64(v_["P-1"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for S = Int64(v_["P+1"]):Int64(v_["NBD"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["1"]):Int64(v_["I-1"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["I+1"]):Int64(v_["NA"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for T = Int64(v_["2"]):Int64(v_["NCE"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            v_["RA"] = v_["RA"*string(Int64(v_["2"]))]
            v_["-RA"] = -1.0*v_["RA"]
            for S = Int64(v_["1"]):Int64(v_["P-1"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for S = Int64(v_["P+1"]):Int64(v_["NBD"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["1"]):Int64(v_["I-1"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["I+1"]):Int64(v_["NA"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for T = Int64(v_["1"]):Int64(v_["1"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            for T = Int64(v_["3"]):Int64(v_["NCE"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            v_["RA"] = v_["RA"*string(Int64(v_["3"]))]
            v_["-RA"] = -1.0*v_["RA"]
            for S = Int64(v_["1"]):Int64(v_["P-1"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for S = Int64(v_["P+1"]):Int64(v_["NBD"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["1"]):Int64(v_["I-1"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["I+1"]):Int64(v_["NA"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for T = Int64(v_["1"]):Int64(v_["2"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            for T = Int64(v_["4"]):Int64(v_["NCE"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            v_["RA"] = v_["RA"*string(Int64(v_["4"]))]
            v_["-RA"] = -1.0*v_["RA"]
            for S = Int64(v_["1"]):Int64(v_["P-1"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for S = Int64(v_["P+1"]):Int64(v_["NBD"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["1"]):Int64(v_["I-1"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["I+1"]):Int64(v_["NA"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for T = Int64(v_["1"]):Int64(v_["3"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            for T = Int64(v_["5"]):Int64(v_["NCE"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            v_["RA"] = v_["RA"*string(Int64(v_["5"]))]
            v_["-RA"] = -1.0*v_["RA"]
            for S = Int64(v_["1"]):Int64(v_["P-1"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for S = Int64(v_["P+1"]):Int64(v_["NBD"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["1"]):Int64(v_["I-1"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["I+1"]):Int64(v_["NA"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for T = Int64(v_["1"]):Int64(v_["4"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            for T = Int64(v_["6"]):Int64(v_["NCE"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            v_["RA"] = v_["RA"*string(Int64(v_["6"]))]
            v_["-RA"] = -1.0*v_["RA"]
            for S = Int64(v_["1"]):Int64(v_["P-1"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for S = Int64(v_["P+1"]):Int64(v_["NBD"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["1"]):Int64(v_["I-1"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["I+1"]):Int64(v_["NA"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B2-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for T = Int64(v_["1"]):Int64(v_["5"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B2-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            v_["P"] = v_["5"]
            v_["P+1"] = 1+v_["P"]
            v_["P-1"] = -1+v_["P"]
            v_["RB"] = v_["RB"*string(I)*","*string(Int64(v_["P"]))]
            v_["-RB"] = -1.0*v_["RB"]
            v_["RA"] = v_["RA"*string(Int64(v_["1"]))]
            v_["-RA"] = -1.0*v_["RA"]
            for S = Int64(v_["1"]):Int64(v_["P-1"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for S = Int64(v_["P+1"]):Int64(v_["NBD"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["1"]):Int64(v_["I-1"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["I+1"]):Int64(v_["NA"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for T = Int64(v_["2"]):Int64(v_["NCE"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            v_["RA"] = v_["RA"*string(Int64(v_["2"]))]
            v_["-RA"] = -1.0*v_["RA"]
            for S = Int64(v_["1"]):Int64(v_["P-1"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for S = Int64(v_["P+1"]):Int64(v_["NBD"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["1"]):Int64(v_["I-1"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["I+1"]):Int64(v_["NA"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for T = Int64(v_["1"]):Int64(v_["1"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            for T = Int64(v_["3"]):Int64(v_["NCE"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            v_["RA"] = v_["RA"*string(Int64(v_["3"]))]
            v_["-RA"] = -1.0*v_["RA"]
            for S = Int64(v_["1"]):Int64(v_["P-1"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for S = Int64(v_["P+1"]):Int64(v_["NBD"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["1"]):Int64(v_["I-1"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["I+1"]):Int64(v_["NA"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for T = Int64(v_["1"]):Int64(v_["2"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            for T = Int64(v_["4"]):Int64(v_["NCE"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            v_["RA"] = v_["RA"*string(Int64(v_["4"]))]
            v_["-RA"] = -1.0*v_["RA"]
            for S = Int64(v_["1"]):Int64(v_["P-1"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for S = Int64(v_["P+1"]):Int64(v_["NBD"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["1"]):Int64(v_["I-1"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["I+1"]):Int64(v_["NA"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for T = Int64(v_["1"]):Int64(v_["3"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            for T = Int64(v_["5"]):Int64(v_["NCE"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            v_["RA"] = v_["RA"*string(Int64(v_["5"]))]
            v_["-RA"] = -1.0*v_["RA"]
            for S = Int64(v_["1"]):Int64(v_["P-1"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for S = Int64(v_["P+1"]):Int64(v_["NBD"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["1"]):Int64(v_["I-1"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["I+1"]):Int64(v_["NA"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for T = Int64(v_["1"]):Int64(v_["4"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            for T = Int64(v_["6"]):Int64(v_["NCE"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            v_["RA"] = v_["RA"*string(Int64(v_["6"]))]
            v_["-RA"] = -1.0*v_["RA"]
            for S = Int64(v_["1"]):Int64(v_["P-1"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for S = Int64(v_["P+1"]):Int64(v_["NBD"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["1"]):Int64(v_["I-1"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["I+1"]):Int64(v_["NA"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for T = Int64(v_["1"]):Int64(v_["5"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            v_["P"] = v_["6"]
            v_["P+1"] = 1+v_["P"]
            v_["P-1"] = -1+v_["P"]
            v_["RB"] = v_["RB"*string(I)*","*string(Int64(v_["P"]))]
            v_["-RB"] = -1.0*v_["RB"]
            v_["RA"] = v_["RA"*string(Int64(v_["1"]))]
            v_["-RA"] = -1.0*v_["RA"]
            for S = Int64(v_["1"]):Int64(v_["P-1"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for S = Int64(v_["P+1"]):Int64(v_["NBD"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["1"]):Int64(v_["I-1"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["I+1"]):Int64(v_["NA"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A1-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for T = Int64(v_["2"]):Int64(v_["NCE"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A1-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["1"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            v_["RA"] = v_["RA"*string(Int64(v_["2"]))]
            v_["-RA"] = -1.0*v_["RA"]
            for S = Int64(v_["1"]):Int64(v_["P-1"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for S = Int64(v_["P+1"]):Int64(v_["NBD"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["1"]):Int64(v_["I-1"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["I+1"]):Int64(v_["NA"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A2-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for T = Int64(v_["1"]):Int64(v_["1"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            for T = Int64(v_["3"]):Int64(v_["NCE"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A2-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["2"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            v_["RA"] = v_["RA"*string(Int64(v_["3"]))]
            v_["-RA"] = -1.0*v_["RA"]
            for S = Int64(v_["1"]):Int64(v_["P-1"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for S = Int64(v_["P+1"]):Int64(v_["NBD"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["1"]):Int64(v_["I-1"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["I+1"]):Int64(v_["NA"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A3-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for T = Int64(v_["1"]):Int64(v_["2"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            for T = Int64(v_["4"]):Int64(v_["NCE"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A3-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["3"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            v_["RA"] = v_["RA"*string(Int64(v_["4"]))]
            v_["-RA"] = -1.0*v_["RA"]
            for S = Int64(v_["1"]):Int64(v_["P-1"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for S = Int64(v_["P+1"]):Int64(v_["NBD"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["1"]):Int64(v_["I-1"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["I+1"]):Int64(v_["NA"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A4-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for T = Int64(v_["1"]):Int64(v_["3"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            for T = Int64(v_["5"]):Int64(v_["NCE"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A4-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["4"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            v_["RA"] = v_["RA"*string(Int64(v_["5"]))]
            v_["-RA"] = -1.0*v_["RA"]
            for S = Int64(v_["1"]):Int64(v_["P-1"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for S = Int64(v_["P+1"]):Int64(v_["NBD"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["1"]):Int64(v_["I-1"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["I+1"]):Int64(v_["NA"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A5-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for T = Int64(v_["1"]):Int64(v_["4"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            for T = Int64(v_["6"]):Int64(v_["NCE"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A5-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["5"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
            v_["RA"] = v_["RA"*string(Int64(v_["6"]))]
            v_["-RA"] = -1.0*v_["RA"]
            for S = Int64(v_["1"]):Int64(v_["P-1"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for S = Int64(v_["P+1"]):Int64(v_["NBD"])
                for R = Int64(v_["1"]):Int64(v_["NA"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["1"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["2"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["3"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["4"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["5"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                    v_["L"] = v_["L"]+v_["1"]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(S)*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["1"]):Int64(v_["I-1"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for R = Int64(v_["I+1"]):Int64(v_["NA"])
                for T = Int64(v_["1"]):Int64(v_["NCE"])
                    v_["L"] = v_["L"]+v_["1"]
                    ig = ig_["I"*string(Int64(v_["L"]))]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A6-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B3-"*string(R)*","*string(Int64(v_["P"]))*","*string(T)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
                end
            end
            for T = Int64(v_["1"]):Int64(v_["5"])
                v_["L"] = v_["L"]+v_["1"]
                ig = ig_["I"*string(Int64(v_["L"]))]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A6-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RA"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(Int64(v_["6"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["RB"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B3-"*string(I)*","*string(Int64(v_["P"]))*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-RB"]))
            end
        end
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+pb.neq+1:pb.m] = zeros(Float64,pb.nge)
        pb.cupper[1:pb.nge] = fill(Inf,pb.nge)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-OOR2-MN-72-1261"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "e_globs"

        pbm = args[1]
        arrset(pbm.efpar,1,0.1e0)
        return pbm

    elseif action == "eA1"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        ALPHA = 0.0e0
        OMEGA = 1.0e0/2.0e0
        OM1 = OMEGA-1.0e0
        OM2 = OMEGA-2.0e0
        CMA = EV_[1]-ALPHA
        BIG = CMA>=pbm.efpar[1]
        if BIG
            F = CMA^OMEGA
        end
        if BIG
            G = OMEGA*CMA^OM1
        end
        if BIG
            H = OMEGA*OM1*CMA^OM2
        end
        if !BIG
            C1 = OMEGA*(2.0e0-OMEGA)*pbm.efpar[1]^OM1
        end
        if !BIG
            C2 = OMEGA*OM1*pbm.efpar[1]^OM2
        end
        if !BIG
            F = ((1.0e0-1.5e0*OMEGA+0.5e0*OMEGA*OMEGA)*pbm.efpar[1]^OMEGA+
             C1*CMA+0.5e0*C2*CMA*CMA)
        end
        if !BIG
            G = C1+C2*CMA
        end
        if !BIG
            H = C2
        end
        f_   = F
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = G
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = H
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eA2"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        ALPHA = 0.0e0
        OMEGA = 2.0e0/3.0e0
        OM1 = OMEGA-1.0e0
        OM2 = OMEGA-2.0e0
        CMA = EV_[1]-ALPHA
        BIG = CMA>=pbm.efpar[1]
        if BIG
            F = CMA^OMEGA
        end
        if BIG
            G = OMEGA*CMA^OM1
        end
        if BIG
            H = OMEGA*OM1*CMA^OM2
        end
        if !BIG
            C1 = OMEGA*(2.0e0-OMEGA)*pbm.efpar[1]^OM1
        end
        if !BIG
            C2 = OMEGA*OM1*pbm.efpar[1]^OM2
        end
        if !BIG
            F = ((1.0e0-1.5e0*OMEGA+0.5e0*OMEGA*OMEGA)*pbm.efpar[1]^OMEGA+
             C1*CMA+0.5e0*C2*CMA*CMA)
        end
        if !BIG
            G = C1+C2*CMA
        end
        if !BIG
            H = C2
        end
        f_   = F
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = G
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = H
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eA3"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        ALPHA = 1.0e0
        OMEGA = 1.0e0/2.0e0
        OM1 = OMEGA-1.0e0
        OM2 = OMEGA-2.0e0
        CMA = EV_[1]-ALPHA
        BIG = CMA>=pbm.efpar[1]
        if BIG
            F = CMA^OMEGA
        end
        if BIG
            G = OMEGA*CMA^OM1
        end
        if BIG
            H = OMEGA*OM1*CMA^OM2
        end
        if !BIG
            C1 = OMEGA*(2.0e0-OMEGA)*pbm.efpar[1]^OM1
        end
        if !BIG
            C2 = OMEGA*OM1*pbm.efpar[1]^OM2
        end
        if !BIG
            F = ((1.0e0-1.5e0*OMEGA+0.5e0*OMEGA*OMEGA)*pbm.efpar[1]^OMEGA+
             C1*CMA+0.5e0*C2*CMA*CMA)
        end
        if !BIG
            G = C1+C2*CMA
        end
        if !BIG
            H = C2
        end
        f_   = F
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = G
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = H
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eA4"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        ALPHA = 1.0e0
        OMEGA = 2.0e0/3.0e0
        OM1 = OMEGA-1.0e0
        OM2 = OMEGA-2.0e0
        CMA = EV_[1]-ALPHA
        BIG = CMA>=pbm.efpar[1]
        if BIG
            F = CMA^OMEGA
        end
        if BIG
            G = OMEGA*CMA^OM1
        end
        if BIG
            H = OMEGA*OM1*CMA^OM2
        end
        if !BIG
            C1 = OMEGA*(2.0e0-OMEGA)*pbm.efpar[1]^OM1
        end
        if !BIG
            C2 = OMEGA*OM1*pbm.efpar[1]^OM2
        end
        if !BIG
            F = ((1.0e0-1.5e0*OMEGA+0.5e0*OMEGA*OMEGA)*pbm.efpar[1]^OMEGA+
             C1*CMA+0.5e0*C2*CMA*CMA)
        end
        if !BIG
            G = C1+C2*CMA
        end
        if !BIG
            H = C2
        end
        f_   = F
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = G
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = H
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eA5"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        ALPHA = 1.5e0
        OMEGA = 1.0e0/2.0e0
        OM1 = OMEGA-1.0e0
        OM2 = OMEGA-2.0e0
        CMA = EV_[1]-ALPHA
        BIG = CMA>=pbm.efpar[1]
        if BIG
            F = CMA^OMEGA
        end
        if BIG
            G = OMEGA*CMA^OM1
        end
        if BIG
            H = OMEGA*OM1*CMA^OM2
        end
        if !BIG
            C1 = OMEGA*(2.0e0-OMEGA)*pbm.efpar[1]^OM1
        end
        if !BIG
            C2 = OMEGA*OM1*pbm.efpar[1]^OM2
        end
        if !BIG
            F = ((1.0e0-1.5e0*OMEGA+0.5e0*OMEGA*OMEGA)*pbm.efpar[1]^OMEGA+
             C1*CMA+0.5e0*C2*CMA*CMA)
        end
        if !BIG
            G = C1+C2*CMA
        end
        if !BIG
            H = C2
        end
        f_   = F
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = G
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = H
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eA6"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        ALPHA = 1.5e0
        OMEGA = 2.0e0/3.0e0
        OM1 = OMEGA-1.0e0
        OM2 = OMEGA-2.0e0
        CMA = EV_[1]-ALPHA
        BIG = CMA>=pbm.efpar[1]
        if BIG
            F = CMA^OMEGA
        end
        if BIG
            G = OMEGA*CMA^OM1
        end
        if BIG
            H = OMEGA*OM1*CMA^OM2
        end
        if !BIG
            C1 = OMEGA*(2.0e0-OMEGA)*pbm.efpar[1]^OM1
        end
        if !BIG
            C2 = OMEGA*OM1*pbm.efpar[1]^OM2
        end
        if !BIG
            F = ((1.0e0-1.5e0*OMEGA+0.5e0*OMEGA*OMEGA)*pbm.efpar[1]^OMEGA+
             C1*CMA+0.5e0*C2*CMA*CMA)
        end
        if !BIG
            G = C1+C2*CMA
        end
        if !BIG
            H = C2
        end
        f_   = F
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = G
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = H
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eB1"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        THETA = 1.0e0/3.0e0
        f_   = EV_[1]^THETA
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = THETA*EV_[1]^(THETA-1.0e0)
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = THETA*(THETA-1.0e0)*EV_[1]^(THETA-2.0e0)
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eB2"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        THETA = 1.0e0/2.0e0
        f_   = EV_[1]^THETA
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = THETA*EV_[1]^(THETA-1.0e0)
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = THETA*(THETA-1.0e0)*EV_[1]^(THETA-2.0e0)
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eB3"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        THETA = 2.0e0/3.0e0
        f_   = EV_[1]^THETA
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = THETA*EV_[1]^(THETA-1.0e0)
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = THETA*(THETA-1.0e0)*EV_[1]^(THETA-2.0e0)
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
            pbm.has_globs = [1,0]
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

