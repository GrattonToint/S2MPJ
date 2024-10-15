function HYDCAR6LS(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HYDCAR6LS
#    *********
# 
#    The hydrocarbon-6 problem by Fletcher.
# 
#    Source: Problem 2a in
#    J.J. More',"A collection of nonlinear model problems"
#    Proceedings of the AMS-SIAM Summer Seminar on the Computational
#    Solution of Nonlinear Systems of Equations, Colorado, 1988.
#    Argonne National Laboratory MCS-P60-0289, 1989.
# 
#    SIF input : N. Gould, Dec 1989
#    Least-squares version of HYDCAR6.SIF, Nick Gould, Jan 2020.
# 
#    classification = "C-SUR2-AN-29-0"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HYDCAR6LS"

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
        v_["N"] = 6
        v_["N-1"] = -1+v_["N"]
        v_["N-2"] = -2+v_["N"]
        v_["M"] = 3
        v_["K"] = 2
        v_["K-"] = -1+v_["K"]
        v_["K+"] = 1+v_["K"]
        v_["A1"] = 9.647
        v_["B1"] = -2998.00
        v_["C1"] = 230.66
        v_["A2"] = 9.953
        v_["B2"] = -3448.10
        v_["C2"] = 235.88
        v_["A3"] = 9.466
        v_["B3"] = -3347.25
        v_["C3"] = 215.31
        v_["AL1"] = 0.0
        v_["ALp1"] = 37.6
        v_["ALpp1"] = 0.0
        v_["AL2"] = 0.0
        v_["ALp2"] = 48.2
        v_["ALpp2"] = 0.0
        v_["AL3"] = 0.0
        v_["ALp3"] = 45.4
        v_["ALpp3"] = 0.0
        v_["BE1"] = 8425.0
        v_["BEp1"] = 24.2
        v_["BEpp1"] = 0.0
        v_["BE2"] = 9395.0
        v_["BEp2"] = 35.6
        v_["BEpp2"] = 0.0
        v_["BE3"] = 10466.0
        v_["BEp3"] = 31.9
        v_["BEpp3"] = 0.0
        v_["FL1"] = 30.0
        v_["FL2"] = 30.0
        v_["FL3"] = 40.0
        v_["FV1"] = 0.0
        v_["FV2"] = 0.0
        v_["FV3"] = 0.0
        v_["TF"] = 100.0
        v_["B"] = 40.0
        v_["D"] = 60.0
        v_["Q"] = 2500000.0
        for I = Int64(v_["0"]):Int64(v_["N-1"])
            v_["PI"*string(I)] = 1.0
        end
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["0"]):Int64(v_["N-1"])
            iv,ix_,_ = s2mpj_ii("T"*string(I),ix_)
            arrset(pb.xnames,iv,"T"*string(I))
            v_["INVPI"*string(I)] = 1.0/v_["PI"*string(I)]
            for J = Int64(v_["1"]):Int64(v_["M"])
                iv,ix_,_ = s2mpj_ii("X"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"X"*string(I)*","*string(J))
            end
        end
        for I = Int64(v_["0"]):Int64(v_["N-2"])
            iv,ix_,_ = s2mpj_ii("V"*string(I),ix_)
            arrset(pb.xnames,iv,"V"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for J = Int64(v_["1"]):Int64(v_["M"])
            ig,ig_,_ = s2mpj_ii("2.1-"*string(J),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(Int64(v_["0"]))*","*string(J)]
            pbm.A[ig,iv] += Float64(v_["B"])
            arrset(pbm.gscale,ig,Float64(1.0e+4))
            ig,ig_,_ = s2mpj_ii("2.3-"*string(J),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(Int64(v_["N-1"]))*","*string(J)]
            pbm.A[ig,iv] += Float64(-1.0)
            for I = Int64(v_["1"]):Int64(v_["N-2"])
                ig,ig_,_ = s2mpj_ii("2.2-"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<>")
                arrset(pbm.gscale,ig,Float64(1.0e+4))
            end
        end
        for I = Int64(v_["0"]):Int64(v_["N-1"])
            ig,ig_,_ = s2mpj_ii("2.7-"*string(I),ig_)
            arrset(gtype,ig,"<>")
        end
        ig,ig_,_ = s2mpj_ii("2.8",ig_)
        arrset(gtype,ig,"<>")
        arrset(pbm.gscale,ig,Float64(1.0e+10))
        for I = Int64(v_["1"]):Int64(v_["N-2"])
            ig,ig_,_ = s2mpj_ii("2.9-"*string(I),ig_)
            arrset(gtype,ig,"<>")
            arrset(pbm.gscale,ig,Float64(1.0e+10))
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        v_["SMALLHF"] = 0.0e+0
        v_["BIGHF"] = 0.0e+0
        for J = Int64(v_["1"]):Int64(v_["M"])
            pbm.gconst[ig_["2.2-"*string(Int64(v_["K"]))*","*string(J)]]  = (
                  Float64(v_["FL"*string(J)]))
            pbm.gconst[ig_["2.2-"*string(Int64(v_["K+"]))*","*string(J)]]  = (
                  Float64(v_["FV"*string(J)]))
            v_["TFTF"] = v_["TF"]*v_["TF"]
            v_["TEMP1"] = v_["TFTF"]*v_["ALpp"*string(J)]
            v_["TEMP2"] = v_["TF"]*v_["ALp"*string(J)]
            v_["TEMP1"] = v_["TEMP1"]+v_["TEMP2"]
            v_["TEMP1"] = v_["TEMP1"]+v_["AL"*string(J)]
            v_["TEMP1"] = v_["TEMP1"]*v_["FL"*string(J)]
            v_["SMALLHF"] = v_["SMALLHF"]+v_["TEMP1"]
            v_["TEMP1"] = v_["TFTF"]*v_["BEpp"*string(J)]
            v_["TEMP2"] = v_["TF"]*v_["BEp"*string(J)]
            v_["TEMP1"] = v_["TEMP1"]+v_["TEMP2"]
            v_["TEMP1"] = v_["TEMP1"]+v_["BE"*string(J)]
            v_["TEMP1"] = v_["TEMP1"]*v_["FV"*string(J)]
            v_["BIGHF"] = v_["BIGHF"]+v_["TEMP1"]
        end
        for I = Int64(v_["0"]):Int64(v_["N-1"])
            pbm.gconst[ig_["2.7-"*string(I)]] = Float64(1.0)
        end
        pbm.gconst[ig_["2.8"]] = Float64(v_["Q"])
        pbm.gconst[ig_["2.9-"*string(Int64(v_["K"]))]] = Float64(v_["SMALLHF"])
        pbm.gconst[ig_["2.9-"*string(Int64(v_["K+"]))]] = Float64(v_["BIGHF"])
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.x0[ix_["X0,1"]] = Float64(0.0)
        pb.x0[ix_["X0,2"]] = Float64(0.2)
        pb.x0[ix_["X0,3"]] = Float64(0.9)
        pb.x0[ix_["X1,1"]] = Float64(0.0)
        pb.x0[ix_["X1,2"]] = Float64(0.2)
        pb.x0[ix_["X1,3"]] = Float64(0.8)
        pb.x0[ix_["X2,1"]] = Float64(0.05)
        pb.x0[ix_["X2,2"]] = Float64(0.3)
        pb.x0[ix_["X2,3"]] = Float64(0.8)
        pb.x0[ix_["X3,1"]] = Float64(0.1)
        pb.x0[ix_["X3,2"]] = Float64(0.3)
        pb.x0[ix_["X3,3"]] = Float64(0.6)
        pb.x0[ix_["X4,1"]] = Float64(0.3)
        pb.x0[ix_["X4,2"]] = Float64(0.5)
        pb.x0[ix_["X4,3"]] = Float64(0.3)
        pb.x0[ix_["X5,1"]] = Float64(0.6)
        pb.x0[ix_["X5,2"]] = Float64(0.6)
        pb.x0[ix_["X5,3"]] = Float64(0.0)
        for I = Int64(v_["0"]):Int64(v_["N-1"])
            pb.x0[ix_["T"*string(I)]] = Float64(100.0)
        end
        for I = Int64(v_["0"]):Int64(v_["N-2"])
            pb.x0[ix_["V"*string(I)]] = Float64(300.0)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "en2PROD", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"P1")
        loaset(elftp,it,2,"P2")
        it,iet_,_ = s2mpj_ii( "ePOLY1PRD", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        loaset(elftp,it,1,"P1")
        loaset(elftp,it,2,"P6")
        loaset(elftp,it,3,"P7")
        loaset(elftp,it,4,"P8")
        it,iet_,_ = s2mpj_ii( "ePOLY2PRD", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        loaset(elftv,it,3,"V3")
        loaset(elftp,it,1,"P1")
        loaset(elftp,it,2,"P2")
        loaset(elftp,it,3,"P6")
        loaset(elftp,it,4,"P7")
        loaset(elftp,it,5,"P8")
        it,iet_,_ = s2mpj_ii( "eEXP2PROD", iet_)
        loaset(elftv,it,1,"V2")
        loaset(elftv,it,2,"V3")
        loaset(elftp,it,1,"P1")
        loaset(elftp,it,2,"P2")
        loaset(elftp,it,3,"P3")
        loaset(elftp,it,4,"P4")
        loaset(elftp,it,5,"P5")
        it,iet_,_ = s2mpj_ii( "eEXP3PROD", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        loaset(elftv,it,3,"V3")
        loaset(elftp,it,1,"P1")
        loaset(elftp,it,2,"P2")
        loaset(elftp,it,3,"P3")
        loaset(elftp,it,4,"P4")
        loaset(elftp,it,5,"P5")
        it,iet_,_ = s2mpj_ii( "eEXP4PROD", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        loaset(elftv,it,3,"V3")
        loaset(elftp,it,1,"P1")
        loaset(elftp,it,2,"P2")
        loaset(elftp,it,3,"P3")
        loaset(elftp,it,4,"P4")
        loaset(elftp,it,5,"P5")
        loaset(elftp,it,6,"P6")
        loaset(elftp,it,7,"P7")
        loaset(elftp,it,8,"P8")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        v_["-D"] = -1.0*v_["D"]
        for J = Int64(v_["1"]):Int64(v_["M"])
            ename = "E11-"*string(J)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"en2PROD")
            arrset(ielftype,ie,iet_["en2PROD"])
            vname = "X"*string(Int64(v_["1"]))*","*string(J)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "V"*string(Int64(v_["0"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(-1.0))
            posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["B"]))
            ename = "E12-"*string(J)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eEXP3PROD")
            arrset(ielftype,ie,iet_["eEXP3PROD"])
            vname = "V"*string(Int64(v_["0"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["0"]))*","*string(J)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "T"*string(Int64(v_["0"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(1.0))
            posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["INVPI"*string(Int64(v_["0"]))]))
            posep = findfirst(x->x=="P3",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["A"*string(J)]))
            posep = findfirst(x->x=="P4",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["B"*string(J)]))
            posep = findfirst(x->x=="P5",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["C"*string(J)]))
            for I = Int64(v_["1"]):Int64(v_["N-2"])
                v_["I-1"] = -1+I
                v_["I+1"] = 1+I
                ename = "E21-"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"en2PROD")
                arrset(ielftype,ie,iet_["en2PROD"])
                vname = "X"*string(Int64(v_["I+1"]))*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "V"*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(-1.0))
                ename = "E22-"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eEXP3PROD")
                arrset(ielftype,ie,iet_["eEXP3PROD"])
                vname = "V"*string(Int64(v_["I-1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(Int64(v_["I-1"]))*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "T"*string(Int64(v_["I-1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(-1.0))
                posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["INVPI"*string(Int64(v_["I-1"]))]))
                posep = findfirst(x->x=="P3",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["A"*string(J)]))
                posep = findfirst(x->x=="P4",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["B"*string(J)]))
                posep = findfirst(x->x=="P5",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["C"*string(J)]))
                ename = "E23-"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"en2PROD")
                arrset(ielftype,ie,iet_["en2PROD"])
                vname = "X"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "V"*string(Int64(v_["I-1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(1.0))
                ename = "E24-"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eEXP3PROD")
                arrset(ielftype,ie,iet_["eEXP3PROD"])
                vname = "V"*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "T"*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(1.0))
                posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["INVPI"*string(I)]))
                posep = findfirst(x->x=="P3",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["A"*string(J)]))
                posep = findfirst(x->x=="P4",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["B"*string(J)]))
                posep = findfirst(x->x=="P5",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["C"*string(J)]))
            end
            for I = Int64(v_["1"]):Int64(v_["K-"])
                ename = "E21-"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["B"]))
                ename = "E23-"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["B"]))
            end
            ename = "E21-"*string(Int64(v_["K"]))*","*string(J)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["-D"]))
            ename = "E23-"*string(Int64(v_["K"]))*","*string(J)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["B"]))
            for I = Int64(v_["K+"]):Int64(v_["N-2"])
                ename = "E21-"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["-D"]))
                ename = "E23-"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["-D"]))
            end
            ename = "E31-"*string(J)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eEXP2PROD")
            arrset(ielftype,ie,iet_["eEXP2PROD"])
            vname = "X"*string(Int64(v_["N-2"]))*","*string(J)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "T"*string(Int64(v_["N-2"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(1.0))
            posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["INVPI"*string(Int64(v_["N-2"]))]))
            posep = findfirst(x->x=="P3",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["A"*string(J)]))
            posep = findfirst(x->x=="P4",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["B"*string(J)]))
            posep = findfirst(x->x=="P5",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["C"*string(J)]))
        end
        for J = Int64(v_["1"]):Int64(v_["M"])
            for I = Int64(v_["0"]):Int64(v_["N-1"])
                ename = "E71-"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eEXP2PROD")
                arrset(ielftype,ie,iet_["eEXP2PROD"])
                vname = "X"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "T"*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(1.0))
                posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["INVPI"*string(I)]))
                posep = findfirst(x->x=="P3",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["A"*string(J)]))
                posep = findfirst(x->x=="P4",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["B"*string(J)]))
                posep = findfirst(x->x=="P5",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["C"*string(J)]))
            end
        end
        for J = Int64(v_["1"]):Int64(v_["M"])
            ename = "E81-"*string(J)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eEXP4PROD")
            arrset(ielftype,ie,iet_["eEXP4PROD"])
            vname = "V"*string(Int64(v_["0"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["0"]))*","*string(J)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "T"*string(Int64(v_["0"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(1.0))
            posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["INVPI"*string(Int64(v_["0"]))]))
            posep = findfirst(x->x=="P3",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["A"*string(J)]))
            posep = findfirst(x->x=="P4",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["B"*string(J)]))
            posep = findfirst(x->x=="P5",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["C"*string(J)]))
            posep = findfirst(x->x=="P6",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["BE"*string(J)]))
            posep = findfirst(x->x=="P7",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["BEp"*string(J)]))
            posep = findfirst(x->x=="P8",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["BEpp"*string(J)]))
            ename = "E82-"*string(J)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ePOLY1PRD")
            arrset(ielftype,ie,iet_["ePOLY1PRD"])
            vname = "X"*string(Int64(v_["0"]))*","*string(J)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "T"*string(Int64(v_["0"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["B"]))
            posep = findfirst(x->x=="P6",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["AL"*string(J)]))
            posep = findfirst(x->x=="P7",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["ALp"*string(J)]))
            posep = findfirst(x->x=="P8",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["ALpp"*string(J)]))
            ename = "E83-"*string(J)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ePOLY2PRD")
            arrset(ielftype,ie,iet_["ePOLY2PRD"])
            vname = "X"*string(Int64(v_["1"]))*","*string(J)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "V"*string(Int64(v_["0"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "T"*string(Int64(v_["1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(-1.0))
            posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["B"]))
            posep = findfirst(x->x=="P6",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["AL"*string(J)]))
            posep = findfirst(x->x=="P7",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["ALp"*string(J)]))
            posep = findfirst(x->x=="P8",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["ALpp"*string(J)]))
            for I = Int64(v_["1"]):Int64(v_["N-2"])
                v_["I-1"] = -1+I
                v_["I+1"] = 1+I
                ename = "E91-"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eEXP4PROD")
                arrset(ielftype,ie,iet_["eEXP4PROD"])
                vname = "V"*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "T"*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(1.0))
                posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["INVPI"*string(I)]))
                posep = findfirst(x->x=="P3",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["A"*string(J)]))
                posep = findfirst(x->x=="P4",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["B"*string(J)]))
                posep = findfirst(x->x=="P5",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["C"*string(J)]))
                posep = findfirst(x->x=="P6",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["BE"*string(J)]))
                posep = findfirst(x->x=="P7",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["BEp"*string(J)]))
                posep = findfirst(x->x=="P8",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["BEpp"*string(J)]))
                ename = "E92-"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"ePOLY2PRD")
                arrset(ielftype,ie,iet_["ePOLY2PRD"])
                vname = "X"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "V"*string(Int64(v_["I-1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "T"*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(1.0))
                posep = findfirst(x->x=="P6",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["AL"*string(J)]))
                posep = findfirst(x->x=="P7",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["ALp"*string(J)]))
                posep = findfirst(x->x=="P8",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["ALpp"*string(J)]))
                ename = "E93-"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eEXP4PROD")
                arrset(ielftype,ie,iet_["eEXP4PROD"])
                vname = "V"*string(Int64(v_["I-1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(Int64(v_["I-1"]))*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "T"*string(Int64(v_["I-1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(-1.0))
                posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["INVPI"*string(Int64(v_["I-1"]))]))
                posep = findfirst(x->x=="P3",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["A"*string(J)]))
                posep = findfirst(x->x=="P4",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["B"*string(J)]))
                posep = findfirst(x->x=="P5",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["C"*string(J)]))
                posep = findfirst(x->x=="P6",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["BE"*string(J)]))
                posep = findfirst(x->x=="P7",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["BEp"*string(J)]))
                posep = findfirst(x->x=="P8",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["BEpp"*string(J)]))
                ename = "E94-"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"ePOLY2PRD")
                arrset(ielftype,ie,iet_["ePOLY2PRD"])
                vname = "X"*string(Int64(v_["I+1"]))*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "V"*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "T"*string(Int64(v_["I+1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(-1.0))
                posep = findfirst(x->x=="P6",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["AL"*string(J)]))
                posep = findfirst(x->x=="P7",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["ALp"*string(J)]))
                posep = findfirst(x->x=="P8",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["ALpp"*string(J)]))
            end
            for I = Int64(v_["1"]):Int64(v_["K-"])
                ename = "E92-"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["B"]))
                ename = "E94-"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["B"]))
            end
            ename = "E92-"*string(Int64(v_["K"]))*","*string(J)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["B"]))
            ename = "E94-"*string(Int64(v_["K"]))*","*string(J)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["-D"]))
            for I = Int64(v_["K+"]):Int64(v_["N-2"])
                ename = "E92-"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["-D"]))
                ename = "E94-"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["-D"]))
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
        for ig in 1:ngrp
            arrset(pbm.grftype,ig,"gL2")
        end
        for J = Int64(v_["1"]):Int64(v_["M"])
            ig = ig_["2.1-"*string(J)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E11-"*string(J)])
            loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["E12-"*string(J)])
            loaset(pbm.grelw,ig,posel, 1.)
            ig = ig_["2.3-"*string(J)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E31-"*string(J)])
            loaset(pbm.grelw,ig,posel,1.)
            ig = ig_["2.8"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E81-"*string(J)])
            loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["E82-"*string(J)])
            loaset(pbm.grelw,ig,posel, 1.)
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E83-"*string(J)])
            loaset(pbm.grelw,ig,posel,1.)
            for I = Int64(v_["1"]):Int64(v_["N-2"])
                ig = ig_["2.2-"*string(I)*","*string(J)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["E21-"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel,1.)
                posel = posel+1
                loaset(pbm.grelt,ig,posel,ie_["E22-"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel, 1.)
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["E23-"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel,1.)
                posel = posel+1
                loaset(pbm.grelt,ig,posel,ie_["E24-"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel, 1.)
                ig = ig_["2.9-"*string(I)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["E91-"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel,1.)
                posel = posel+1
                loaset(pbm.grelt,ig,posel,ie_["E92-"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel, 1.)
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["E93-"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel,1.)
                posel = posel+1
                loaset(pbm.grelt,ig,posel,ie_["E94-"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel, 1.)
            end
            for I = Int64(v_["0"]):Int64(v_["N-1"])
                ig = ig_["2.7-"*string(I)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["E71-"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel,1.)
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        pb.objlower = 0.0
#    Solution
# LO SOLTN               0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-SUR2-AN-29-0"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
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
        f_   = pbm.elpar[iel_][1]*EV_[1]*(EV_[2]+pbm.elpar[iel_][2])
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = pbm.elpar[iel_][1]*(EV_[2]+pbm.elpar[iel_][2])
            g_[2] = pbm.elpar[iel_][1]*EV_[1]
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = pbm.elpar[iel_][1]
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

    elseif action == "ePOLY1PRD"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        POLY  = (
              pbm.elpar[iel_][2]+pbm.elpar[iel_][3]*EV_[2]+pbm.elpar[iel_][4]*EV_[2]*EV_[2])
        DPOLY = pbm.elpar[iel_][3]+2.0*pbm.elpar[iel_][4]*EV_[2]
        f_   = pbm.elpar[iel_][1]*EV_[1]*POLY
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = pbm.elpar[iel_][1]*POLY
            g_[2] = pbm.elpar[iel_][1]*EV_[1]*DPOLY
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = pbm.elpar[iel_][1]*DPOLY
                H_[2,1] = H_[1,2]
                H_[2,2] = pbm.elpar[iel_][1]*EV_[1]*2.0e+0*pbm.elpar[iel_][4]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "ePOLY2PRD"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        POLY  = (
              pbm.elpar[iel_][3]+pbm.elpar[iel_][4]*EV_[3]+pbm.elpar[iel_][5]*EV_[3]*EV_[3])
        DPOLY = pbm.elpar[iel_][4]+2.0*pbm.elpar[iel_][5]*EV_[3]
        f_   = pbm.elpar[iel_][1]*EV_[1]*(pbm.elpar[iel_][2]+EV_[2])*POLY
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = pbm.elpar[iel_][1]*(pbm.elpar[iel_][2]+EV_[2])*POLY
            g_[2] = pbm.elpar[iel_][1]*EV_[1]*POLY
            g_[3] = pbm.elpar[iel_][1]*EV_[1]*(pbm.elpar[iel_][2]+EV_[2])*DPOLY
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,2] = pbm.elpar[iel_][1]*POLY
                H_[2,1] = H_[1,2]
                H_[1,3] = pbm.elpar[iel_][1]*(pbm.elpar[iel_][2]+EV_[2])*DPOLY
                H_[3,1] = H_[1,3]
                H_[2,3] = pbm.elpar[iel_][1]*EV_[1]*DPOLY
                H_[3,2] = H_[2,3]
                H_[3,3]  = (
                      pbm.elpar[iel_][1]*EV_[1]*(pbm.elpar[iel_][2]+EV_[2])*2.0e+0*pbm.elpar[iel_][5])
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eEXP2PROD"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        EXPROD  = (
              pbm.elpar[iel_][1]*pbm.elpar[iel_][2]*exp(pbm.elpar[iel_][3]+(pbm.elpar[iel_][4]/(EV_[2]+pbm.elpar[iel_][5]))))
        F = EV_[1]*EXPROD
        f_   = F
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EXPROD
            g_[2] = -EV_[1]*EXPROD*pbm.elpar[iel_][4]/(pbm.elpar[iel_][5]+EV_[2])^2
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = -EXPROD*pbm.elpar[iel_][4]/(pbm.elpar[iel_][5]+EV_[2])^2
                H_[2,1] = H_[1,2]
                H_[2,2] = (F*(pbm.elpar[iel_][4]/(pbm.elpar[iel_][5]+EV_[2])^2)^2+
                     2.0e+0*F*pbm.elpar[iel_][4]/(pbm.elpar[iel_][5]+EV_[2])^3)
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eEXP3PROD"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        EXPROD  = (
              pbm.elpar[iel_][1]*pbm.elpar[iel_][2]*exp(pbm.elpar[iel_][3]+(pbm.elpar[iel_][4]/(EV_[3]+pbm.elpar[iel_][5]))))
        F = EV_[1]*EV_[2]*EXPROD
        TERM = -pbm.elpar[iel_][4]/(pbm.elpar[iel_][5]+EV_[3])^2
        f_   = F
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[2]*EXPROD
            g_[2] = EV_[1]*EXPROD
            g_[3] = F*TERM
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,2] = EXPROD
                H_[2,1] = H_[1,2]
                H_[1,3] = EV_[2]*EXPROD*TERM
                H_[3,1] = H_[1,3]
                H_[2,3] = EV_[1]*EXPROD*TERM
                H_[3,2] = H_[2,3]
                H_[3,3]  = (
                      F*(TERM*TERM+2.0e+0*pbm.elpar[iel_][4]/(pbm.elpar[iel_][5]+EV_[3])^3))
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eEXP4PROD"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        EXPROD  = (
              pbm.elpar[iel_][1]*pbm.elpar[iel_][2]*exp(pbm.elpar[iel_][3]+(pbm.elpar[iel_][4]/(EV_[3]+pbm.elpar[iel_][5]))))
        F = EV_[1]*EV_[2]*EXPROD
        POLY  = (
              pbm.elpar[iel_][6]+pbm.elpar[iel_][7]*EV_[3]+pbm.elpar[iel_][8]*EV_[3]*EV_[3])
        DPOLY = pbm.elpar[iel_][7]+2.0*pbm.elpar[iel_][8]*EV_[3]
        TERM = DPOLY-POLY*pbm.elpar[iel_][4]/(pbm.elpar[iel_][5]+EV_[3])^2
        f_   = F*POLY
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[2]*EXPROD*POLY
            g_[2] = EV_[1]*EXPROD*POLY
            g_[3] = F*TERM
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,2] = EXPROD*POLY
                H_[2,1] = H_[1,2]
                H_[1,3] = EV_[2]*EXPROD*TERM
                H_[3,1] = H_[1,3]
                H_[2,3] = EV_[1]*EXPROD*TERM
                H_[3,2] = H_[2,3]
                H_[3,3]  = (
                      F*(-(pbm.elpar[iel_][4]/(pbm.elpar[iel_][5]+EV_[3])^2)*TERM+2.0*pbm.elpar[iel_][8]-DPOLY*pbm.elpar[iel_][4]/(pbm.elpar[iel_][5]+EV_[3])^2+2.0e+0*POLY*pbm.elpar[iel_][4]/(pbm.elpar[iel_][5]+EV_[3])^3))
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

