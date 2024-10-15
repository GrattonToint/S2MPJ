function MRIBASIS(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ******** 
# 
#    An optmization problem arising in the design of medical apparatus.
# 
#    Source:
#    Contribution from a LANCELOT user.
# 
#    SIF input: Arie Quist, TU Delft (NL), 1994.
#    Adaptation for CUTE: Ph. Toint, November 1994.
# 
#    classification = "C-LOR2-MY-36-55"
# 
#    useful constants
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "MRIBASIS"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["0"] = 0
        v_["1"] = 1
        v_["2"] = 2
        v_["3"] = 3
        v_["R2"] = Float64(v_["2"])
        v_["k1"] = 22.8443
        v_["k2"] = 12.4402
        v_["k3"] = 5.23792
        v_["k4"] = 5.12238
        v_["k5"] = 6.44999
        v_["k6"] = 5.32383
        v_["k7"] = 0.58392
        v_["k8"] = 3.94584
        v_["k9"] = -2.75584
        v_["k10"] = 32.0669
        v_["k11"] = 18.2179
        v_["k12"] = 31.7496
        v_["S1,1"] = -0.377126
        v_["S1,2"] = 0.919679
        v_["S1,3"] = 0.109389
        v_["S2,1"] = 0.634857
        v_["S2,2"] = 0.170704
        v_["S2,3"] = 0.753536
        v_["S3,1"] = 0.674338
        v_["S3,2"] = 0.353624
        v_["S3,3"] = -0.648242
        v_["xlo1"] = 1*v_["k7"]
        v_["k10/2"] = v_["k10"]/v_["R2"]
        v_["k8/2"] = v_["k8"]/v_["R2"]
        v_["xup1"] = v_["k10/2"]-v_["k8/2"]
        v_["xlo2"] = v_["k10/2"]+v_["k8/2"]
        v_["k13"] = v_["k7"]*v_["k5"]
        v_["k14"] = v_["k8/2"]*v_["k6"]
        v_["xup2"] = 1*v_["k12"]
        v_["-k1"] = -1*v_["k1"]
        v_["-k2"] = -1*v_["k2"]
        v_["k4-"] = v_["k4"]-v_["k14"]
        v_["-k3"] = -1*v_["k3"]
        v_["-S1,1"] = -1*v_["S1,1"]
        v_["-S1,2"] = -1*v_["S1,2"]
        v_["-S1,3"] = -1*v_["S1,3"]
        v_["-S2,1"] = -1*v_["S2,1"]
        v_["-S2,2"] = -1*v_["S2,2"]
        v_["-S2,3"] = -1*v_["S2,3"]
        v_["-S3,1"] = -1*v_["S3,1"]
        v_["-S3,2"] = -1*v_["S3,2"]
        v_["-S3,3"] = -1*v_["S3,3"]
        v_["2S1,1"] = 2*v_["S1,1"]
        v_["2S1,2"] = 2*v_["S1,2"]
        v_["2S1,3"] = 2*v_["S1,3"]
        v_["2S2,1"] = 2*v_["S2,1"]
        v_["2S2,2"] = 2*v_["S2,2"]
        v_["2S2,3"] = 2*v_["S2,3"]
        v_["2S3,1"] = 2*v_["S3,1"]
        v_["2S3,2"] = 2*v_["S3,2"]
        v_["2S3,3"] = 2*v_["S3,3"]
        v_["-2S1,1"] = -2*v_["S1,1"]
        v_["-2S1,2"] = -2*v_["S1,2"]
        v_["-2S1,3"] = -2*v_["S1,3"]
        v_["-2S2,1"] = -2*v_["S2,1"]
        v_["-2S2,2"] = -2*v_["S2,2"]
        v_["-2S2,3"] = -2*v_["S2,3"]
        v_["-2S3,1"] = -2*v_["S3,1"]
        v_["-2S3,2"] = -2*v_["S3,2"]
        v_["-2S3,3"] = -2*v_["S3,3"]
        v_["Llo1,1"] = v_["S1,1"]*v_["k5"]
        v_["Llo1,2"] = v_["S1,2"]*v_["k5"]
        v_["Llo1,3"] = v_["S1,3"]*v_["k5"]
        v_["Lup1,1"] = v_["S1,1"]*v_["k6"]
        v_["Lup1,2"] = v_["S1,2"]*v_["k6"]
        v_["Lup1,3"] = v_["S1,3"]*v_["k6"]
        v_["Llo2,1"] = v_["S1,1"]*v_["k6"]
        v_["Llo2,2"] = v_["S1,2"]*v_["k6"]
        v_["Llo2,3"] = v_["S1,3"]*v_["k6"]
        v_["4"] = 4
        v_["xm"] = 6
        v_["Lm"] = 4
        v_["xm-"] = -1+v_["xm"]
        v_["xm-2"] = -2+v_["xm"]
        v_["Lm-"] = -1+v_["Lm"]
        v_["Lm-2"] = -2+v_["Lm"]
        v_["R12"] = 12
        v_["1/12"] = 1/v_["R12"]
        v_["-1/12"] = -1*v_["1/12"]
        v_["R0"] = Float64(v_["0"])
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for j = Int64(v_["1"]):Int64(v_["2"])
            for k = Int64(v_["1"]):Int64(v_["xm"])
                iv,ix_,_ = s2mpj_ii("X"*string(j)*","*string(k),ix_)
                arrset(pb.xnames,iv,"X"*string(j)*","*string(k))
            end
        end
        for i = Int64(v_["1"]):Int64(v_["3"])
            for j = Int64(v_["1"]):Int64(v_["2"])
                for k = Int64(v_["1"]):Int64(v_["Lm"])
                    iv,ix_,_ = s2mpj_ii("L"*string(i)*","*string(j)*","*string(k),ix_)
                    arrset(pb.xnames,iv,"L"*string(i)*","*string(j)*","*string(k))
                end
            end
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("Object",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X"*string(Int64(v_["2"]))*","*string(Int64(v_["xm"]))]
        pbm.A[ig,iv] += Float64(1)
        for j = Int64(v_["1"]):Int64(v_["2"])
            for k = Int64(v_["2"]):Int64(v_["xm-2"])
                v_["k+"] = 1+k
                ig,ig_,_ = s2mpj_ii("PS"*string(j)*","*string(k),ig_)
                arrset(gtype,ig,">=")
                arrset(pb.cnames,ig,"PS"*string(j)*","*string(k))
                iv = ix_["X"*string(j)*","*string(Int64(v_["k+"]))]
                pbm.A[ig,iv] += Float64(1)
                iv = ix_["X"*string(j)*","*string(k)]
                pbm.A[ig,iv] += Float64(-1)
            end
        end
        ig,ig_,_ = s2mpj_ii("PL",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"PL")
        iv = ix_["X"*string(Int64(v_["2"]))*","*string(Int64(v_["xm"]))]
        pbm.A[ig,iv] += Float64(1)
        iv = ix_["X"*string(Int64(v_["2"]))*","*string(Int64(v_["xm-"]))]
        pbm.A[ig,iv] += Float64(-1)
        for i = Int64(v_["1"]):Int64(v_["3"])
            for j = Int64(v_["1"]):Int64(v_["2"])
                for k = Int64(v_["1"]):Int64(v_["Lm-"])
                    v_["k+"] = 1+k
                    v_["2k"] = 2*k
                    v_["2k-"] = -1+v_["2k"]
                    ig,ig_,_ = s2mpj_ii("SU"*string(i)*","*string(j)*","*string(k),ig_)
                    arrset(gtype,ig,"<=")
                    arrset(pb.cnames,ig,"SU"*string(i)*","*string(j)*","*string(k))
                    iv = ix_["L"*string(i)*","*string(j)*","*string(Int64(v_["k+"]))]
                    pbm.A[ig,iv] += Float64(1)
                    iv = ix_["L"*string(i)*","*string(j)*","*string(k)]
                    pbm.A[ig,iv] += Float64(-1)
                    iv = ix_["X"*string(j)*","*string(Int64(v_["2k"]))]
                    pbm.A[ig,iv] += Float64(v_["-k1"])
                    iv = ix_["X"*string(j)*","*string(Int64(v_["2k-"]))]
                    pbm.A[ig,iv] += Float64(v_["k1"])
                    ig,ig_,_ = s2mpj_ii("SL"*string(i)*","*string(j)*","*string(k),ig_)
                    arrset(gtype,ig,">=")
                    arrset(pb.cnames,ig,"SL"*string(i)*","*string(j)*","*string(k))
                    iv = ix_["L"*string(i)*","*string(j)*","*string(Int64(v_["k+"]))]
                    pbm.A[ig,iv] += Float64(1)
                    iv = ix_["L"*string(i)*","*string(j)*","*string(k)]
                    pbm.A[ig,iv] += Float64(-1)
                    iv = ix_["X"*string(j)*","*string(Int64(v_["2k"]))]
                    pbm.A[ig,iv] += Float64(v_["k1"])
                    iv = ix_["X"*string(j)*","*string(Int64(v_["2k-"]))]
                    pbm.A[ig,iv] += Float64(v_["-k1"])
                end
            end
        end
        for i = Int64(v_["1"]):Int64(v_["3"])
            for k = Int64(v_["2"]):Int64(v_["Lm-"])
                ig,ig_,_ = s2mpj_ii("cc1",ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"cc1")
                iv = ix_["L"*string(i)*","*string(Int64(v_["2"]))*","*string(k)]
                pbm.A[ig,iv] += Float64(v_["S"*string(Int64(v_["3"]))*","*string(i)])
            end
        end
        for i = Int64(v_["1"]):Int64(v_["3"])
            ig,ig_,_ = s2mpj_ii("c2const"*string(i),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"c2const"*string(i))
            iv  = (
                  ix_["L"*string(i)*","*string(Int64(v_["2"]))*","*string(Int64(v_["Lm"]))])
            pbm.A[ig,iv] += Float64(v_["k10"])
        end
        ig,ig_,_ = s2mpj_ii("c3con1",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"c3con1")
        ig,ig_,_ = s2mpj_ii("c3con2",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"c3con2")
        ig,ig_,_ = s2mpj_ii("c4const",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"c4const")
        for i = Int64(v_["1"]):Int64(v_["3"])
            ig,ig_,_ = s2mpj_ii("c5con"*string(i),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"c5con"*string(i))
        end
        for i = Int64(v_["1"]):Int64(v_["2"])
            ig,ig_,_ = s2mpj_ii("c6cn"*string(i),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"c6cn"*string(i))
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
        v_["Opmr1"] = v_["k3"]*v_["S2,1"]
        v_["Opmr2"] = v_["k3"]*v_["S2,2"]
        v_["Opmr3"] = v_["k3"]*v_["S2,3"]
        for i = Int64(v_["1"]):Int64(v_["3"])
            pbm.gconst[ig_["c2const"*string(i)]] = Float64(v_["Opmr"*string(i)])
        end
        pbm.gconst[ig_["c3con1"]] = Float64(v_["k4-"])
        pbm.gconst[ig_["c3con2"]] = Float64(v_["k13"])
        pbm.gconst[ig_["c5con1"]] = Float64(v_["k13"])
        pbm.gconst[ig_["c5con2"]] = Float64(v_["-k3"])
        pbm.gconst[ig_["c5con3"]] = Float64(v_["k9"])
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        for j = Int64(v_["1"]):Int64(v_["2"])
            for k = Int64(v_["2"]):Int64(v_["xm-"])
                pb.xlower[ix_["X"*string(j)*","*string(k)]] = v_["xlo"*string(j)]
                pb.xupper[ix_["X"*string(j)*","*string(k)]] = v_["xup"*string(j)]
            end
            pb.xlower[ix_["X"*string(j)*","*string(Int64(v_["1"]))]]  = (
                  v_["xlo"*string(j)])
            pb.xupper[ix_["X"*string(j)*","*string(Int64(v_["1"]))]]  = (
                  v_["xlo"*string(j)])
        end
        pb.xlower[ix_["X"*string(Int64(v_["1"]))*","*string(Int64(v_["xm"]))]]  = (
              v_["xup1"])
        pb.xupper[ix_["X"*string(Int64(v_["1"]))*","*string(Int64(v_["xm"]))]]  = (
              v_["xup1"])
        pb.xlower[ix_["X"*string(Int64(v_["2"]))*","*string(Int64(v_["xm"]))]]  = (
              v_["k11"])
        pb.xupper[ix_["X"*string(Int64(v_["2"]))*","*string(Int64(v_["xm"]))]]  = (
              v_["k12"])
        for i = Int64(v_["1"]):Int64(v_["3"])
            for j = Int64(v_["1"]):Int64(v_["2"])
                for k = Int64(v_["2"]):Int64(v_["Lm-"])
                    pb.xlower[ix_["L"*string(i)*","*string(j)*","*string(k)]] = v_["-k2"]
                    pb.xupper[ix_["L"*string(i)*","*string(j)*","*string(k)]] = v_["k2"]
                end
                pb.xlower[ix_["L"*string(i)*","*string(j)*","*string(Int64(v_["1"]))]]  = (
                      v_["Llo"*string(j)*","*string(i)])
                pb.xupper[ix_["L"*string(i)*","*string(j)*","*string(Int64(v_["1"]))]]  = (
                      v_["Llo"*string(j)*","*string(i)])
            end
            pb.xlower[ix_["L"*string(i)*","*string(Int64(v_["1"]))*","*string(Int64(v_["Lm"]))]] = v_["Lup"*string(Int64(v_["1"]))*","*string(i)]
            pb.xupper[ix_["L"*string(i)*","*string(Int64(v_["1"]))*","*string(Int64(v_["Lm"]))]] = v_["Lup"*string(Int64(v_["1"]))*","*string(i)]
            pb.xlower[ix_["L"*string(i)*","*string(Int64(v_["2"]))*","*string(Int64(v_["Lm"]))]] = v_["-k2"]
            pb.xupper[ix_["L"*string(i)*","*string(Int64(v_["2"]))*","*string(Int64(v_["Lm"]))]] = v_["k2"]
        end
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        v_["intlen1"] = v_["xup1"]-v_["xlo1"]
        v_["Rxm-"] = Float64(v_["xm-"])
        v_["dx1"] = v_["intlen1"]/v_["Rxm-"]
        v_["intlen2"] = v_["k11"]-v_["xlo2"]
        v_["dx2"] = v_["intlen2"]/v_["Rxm-"]
        for k = Int64(v_["1"]):Int64(v_["xm-"])
            v_["Rk"] = Float64(k)
            v_["dist1"] = v_["dx1"]*v_["Rk"]
            v_["strtv1"] = v_["xlo1"]+v_["dist1"]
            v_["dist2"] = v_["dx2"]*v_["Rk"]
            v_["strtv2"] = v_["xlo2"]+v_["dist2"]
            v_["k+"] = 1+k
            pb.x0[ix_["X"*string(Int64(v_["1"]))*","*string(Int64(v_["k+"]))]]  = (
                  Float64(v_["strtv1"]))
            pb.x0[ix_["X"*string(Int64(v_["2"]))*","*string(Int64(v_["k+"]))]]  = (
                  Float64(v_["strtv2"]))
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "euv1", iet_)
        loaset(elftv,it,1,"v1")
        loaset(elftv,it,2,"v2")
        loaset(elftv,it,3,"v3")
        loaset(elftv,it,4,"v4")
        it,iet_,_ = s2mpj_ii( "euv2", iet_)
        loaset(elftv,it,1,"v1")
        loaset(elftv,it,2,"v2")
        loaset(elftv,it,3,"v3")
        it,iet_,_ = s2mpj_ii( "euvw1", iet_)
        loaset(elftv,it,1,"v1")
        loaset(elftv,it,2,"v2")
        loaset(elftv,it,3,"v3")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"p1")
        it,iet_,_ = s2mpj_ii( "emo", iet_)
        loaset(elftv,it,1,"s1")
        loaset(elftv,it,2,"s2")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for i = Int64(v_["1"]):Int64(v_["3"])
            for j = Int64(v_["1"]):Int64(v_["2"])
                for k = Int64(v_["1"]):Int64(v_["Lm-"])
                    v_["2k"] = 2*k
                    v_["k+"] = 1+k
                    v_["2k-"] = -1+v_["2k"]
                    ename = "e1"*string(i)*","*string(j)*","*string(k)
                    ie,ie_,_  = s2mpj_ii(ename,ie_)
                    arrset(pbm.elftype,ie,"euv1")
                    arrset(ielftype,ie,iet_["euv1"])
                    vname = "X"*string(j)*","*string(Int64(v_["2k"]))
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="v1",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    vname = "X"*string(j)*","*string(Int64(v_["2k-"]))
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="v2",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    vname = "L"*string(i)*","*string(j)*","*string(Int64(v_["k+"]))
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="v3",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    vname = "L"*string(i)*","*string(j)*","*string(k)
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="v4",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    ename = "e3"*string(i)*","*string(j)*","*string(k)
                    ie,ie_,_  = s2mpj_ii(ename,ie_)
                    arrset(pbm.elftype,ie,"euvw1")
                    arrset(ielftype,ie,iet_["euvw1"])
                    vname = "L"*string(i)*","*string(j)*","*string(k)
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="v1",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    vname = "X"*string(j)*","*string(Int64(v_["2k"]))
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="v2",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    vname = "X"*string(j)*","*string(Int64(v_["2k-"]))
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="v3",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    posep = findfirst(x->x=="p1",elftp[ielftype[ie]])
                    loaset(pbm.elpar,ie,posep,Float64(v_["1/12"]))
                    ename = "e5"*string(i)*","*string(j)*","*string(k)
                    ie,ie_,_  = s2mpj_ii(ename,ie_)
                    arrset(pbm.elftype,ie,"euvw1")
                    arrset(ielftype,ie,iet_["euvw1"])
                    vname = "L"*string(i)*","*string(j)*","*string(Int64(v_["k+"]))
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="v1",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    vname = "X"*string(j)*","*string(Int64(v_["2k"]))
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="v2",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    vname = "X"*string(j)*","*string(Int64(v_["2k-"]))
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="v3",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    posep = findfirst(x->x=="p1",elftp[ielftype[ie]])
                    loaset(pbm.elpar,ie,posep,Float64(v_["-1/12"]))
                end
                for k = Int64(v_["1"]):Int64(v_["Lm-2"])
                    v_["2k"] = 2*k
                    v_["k+"] = 1+k
                    v_["2k+"] = 1+v_["2k"]
                    ename = "e2"*string(i)*","*string(j)*","*string(k)
                    ie,ie_,_  = s2mpj_ii(ename,ie_)
                    arrset(pbm.elftype,ie,"euv2")
                    arrset(ielftype,ie,iet_["euv2"])
                    vname = "X"*string(j)*","*string(Int64(v_["2k+"]))
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="v1",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    vname = "X"*string(j)*","*string(Int64(v_["2k"]))
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="v2",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    vname = "L"*string(i)*","*string(j)*","*string(Int64(v_["k+"]))
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="v3",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    ename = "e4"*string(i)*","*string(j)*","*string(k)
                    ie,ie_,_  = s2mpj_ii(ename,ie_)
                    arrset(pbm.elftype,ie,"euvw1")
                    arrset(ielftype,ie,iet_["euvw1"])
                    vname = "L"*string(i)*","*string(j)*","*string(Int64(v_["k+"]))
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="v1",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    vname = "X"*string(j)*","*string(Int64(v_["2k+"]))
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="v2",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    vname = "X"*string(j)*","*string(Int64(v_["2k"]))
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="v3",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    posep = findfirst(x->x=="p1",elftp[ielftype[ie]])
                    loaset(pbm.elpar,ie,posep,Float64(v_["R0"]))
                end
            end
        end
        for i = Int64(v_["1"]):Int64(v_["3"])
            ename = "factr"*string(i)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"emo")
            arrset(ielftype,ie,iet_["emo"])
            vname = "X"*string(Int64(v_["2"]))*","*string(Int64(v_["xm"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="s1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "L"*string(i)*","*string(Int64(v_["2"]))*","*string(Int64(v_["Lm"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="s2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for i = Int64(v_["1"]):Int64(v_["3"])
            ig = ig_["c2const"*string(i)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["factr"*string(i)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        for j = Int64(v_["1"]):Int64(v_["2"])
            for i = Int64(v_["1"]):Int64(v_["3"])
                for k = Int64(v_["1"]):Int64(v_["Lm-"])
                    ig = ig_["c3con"*string(j)]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["e1"*string(i)*","*string(Int64(v_["2"]))*","*string(k)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["S"*string(Int64(v_["1"]))*","*string(i)]))
                end
                for k = Int64(v_["1"]):Int64(v_["Lm-2"])
                    ig = ig_["c3con"*string(j)]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["e2"*string(i)*","*string(Int64(v_["2"]))*","*string(k)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["S"*string(Int64(v_["1"]))*","*string(i)]))
                end
            end
        end
        for i = Int64(v_["1"]):Int64(v_["3"])
            for k = Int64(v_["1"]):Int64(v_["Lm-"])
                ig = ig_["c4const"]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["e1"*string(i)*","*string(Int64(v_["2"]))*","*string(k)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["S"*string(Int64(v_["2"]))*","*string(i)]))
            end
            for k = Int64(v_["1"]):Int64(v_["Lm-2"])
                ig = ig_["c4const"]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["e2"*string(i)*","*string(Int64(v_["2"]))*","*string(k)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["S"*string(Int64(v_["2"]))*","*string(i)]))
            end
        end
        for j = Int64(v_["1"]):Int64(v_["3"])
            for i = Int64(v_["1"]):Int64(v_["3"])
                for k = Int64(v_["1"]):Int64(v_["Lm-"])
                    ig = ig_["c5con"*string(j)]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["e1"*string(i)*","*string(Int64(v_["1"]))*","*string(k)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-S"*string(j)*","*string(i)]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["e1"*string(i)*","*string(Int64(v_["2"]))*","*string(k)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["S"*string(j)*","*string(i)]))
                end
                for k = Int64(v_["1"]):Int64(v_["Lm-2"])
                    ig = ig_["c5con"*string(j)]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["e2"*string(i)*","*string(Int64(v_["1"]))*","*string(k)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-S"*string(j)*","*string(i)]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["e2"*string(i)*","*string(Int64(v_["2"]))*","*string(k)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["S"*string(j)*","*string(i)]))
                end
            end
        end
        for j = Int64(v_["1"]):Int64(v_["2"])
            for i = Int64(v_["1"]):Int64(v_["3"])
                for k = Int64(v_["1"]):Int64(v_["Lm-"])
                    ig = ig_["c6cn"*string(j)]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["e3"*string(i)*","*string(Int64(v_["1"]))*","*string(k)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["S"*string(j)*","*string(i)]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["e5"*string(i)*","*string(Int64(v_["1"]))*","*string(k)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["S"*string(j)*","*string(i)]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["e3"*string(i)*","*string(Int64(v_["2"]))*","*string(k)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-S"*string(j)*","*string(i)]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["e5"*string(i)*","*string(Int64(v_["2"]))*","*string(k)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-S"*string(j)*","*string(i)]))
                end
                for k = Int64(v_["1"]):Int64(v_["Lm-2"])
                    ig = ig_["c6cn"*string(j)]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["e4"*string(i)*","*string(Int64(v_["1"]))*","*string(k)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["2S"*string(j)*","*string(i)]))
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["e4"*string(i)*","*string(Int64(v_["2"]))*","*string(k)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,Float64(v_["-2S"*string(j)*","*string(i)]))
                end
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               18.2179000000
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.clower[pb.nle+pb.neq+1:pb.m] = zeros(Float64,pb.nge)
        pb.cupper[1:pb.nge] = fill(Inf,pb.nge)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-LOR2-MY-36-55"
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

    elseif action == "euv1"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,2,4)
        IV_ =  zeros(Float64,2)
        U_[1,1] = U_[1,1]+1
        U_[1,2] = U_[1,2]-1
        U_[2,3] = U_[2,3]+1
        U_[2,4] = U_[2,4]+1
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        f_   = 0.5e0*IV_[1]*IV_[2]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 0.5e0*IV_[2]
            g_[2] = 0.5e0*IV_[1]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = 0.5e0
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

    elseif action == "euv2"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,2,3)
        IV_ =  zeros(Float64,2)
        U_[1,1] = U_[1,1]+1
        U_[1,2] = U_[1,2]-1
        U_[2,3] = U_[2,3]+1
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
                H_[1,2] = 1.0e0
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

    elseif action == "euvw1"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        d14 = 0.25e0
        f_    = (
              EV_[1]*(EV_[2]-EV_[3])*((d14+pbm.elpar[iel_][1])*EV_[2]+(d14-pbm.elpar[iel_][1])*EV_[3]))
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1]  = (
                  (EV_[2]-EV_[3])*((d14+pbm.elpar[iel_][1])*EV_[2]+(d14-pbm.elpar[iel_][1])*EV_[3]))
            g_[2]  = (
                  EV_[1]*2.0e0*(EV_[2]*(d14+pbm.elpar[iel_][1])-pbm.elpar[iel_][1]*EV_[3]))
            g_[3]  = (
                  EV_[1]*2.0e0*(-EV_[3]*(d14-pbm.elpar[iel_][1])-pbm.elpar[iel_][1]*EV_[2]))
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,2] = 2.0e0*(EV_[2]*(d14+pbm.elpar[iel_][1])-pbm.elpar[iel_][1]*EV_[3])
                H_[2,1] = H_[1,2]
                H_[1,3] = 2.0e0*(-EV_[3]*(d14-pbm.elpar[iel_][1])-pbm.elpar[iel_][1]*EV_[2])
                H_[3,1] = H_[1,3]
                H_[2,2] = EV_[1]*2.0e0*(d14+pbm.elpar[iel_][1])
                H_[2,3] = -EV_[1]*2.0e0*pbm.elpar[iel_][1]
                H_[3,2] = H_[2,3]
                H_[3,3] = -EV_[1]*2.0e0*(d14-pbm.elpar[iel_][1])
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "emo"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = -EV_[2]*EV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = -EV_[2]
            g_[2] = -EV_[1]
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = -1.0e0
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

