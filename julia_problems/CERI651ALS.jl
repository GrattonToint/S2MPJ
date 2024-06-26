function CERI651ALS(action,args...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CERI651ALS
#    *********
# 
#    ISIS Data fitting problem CERI651A given as an inconsistent set of
#    nonlinear equations.
# 
#    Fit: y = c + l * x + I*A*B/2(A+B) *
#               [ exp( A*[A*S^2+2(x-X0)]/2) * erfc( A*S^2+(x-X0)/S*sqrt(2) ) +
#                 exp( B*[B*S^2+2(x-X0)]/2) * erfc( B*S^2+(x-X0)/S*sqrt(2) ) ]
# 
#    Source: fit to a sum of a linear background and a back-to-back exponential
#    using data enginx_ceria193749_spectrum_number_651_vana_corrected-0
#    from Mantid (http://www.mantidproject.org)
# 
#    subset X in [36844.7449265, 37300.5256846]
# 
#    SIF input: Nick Gould and Tyrone Rees, Mar 2016
#    Least-squares version of CERI651A.SIF, Nick Gould, Jan 2020.
# 
#    classification = "SUR2-MN-7-0"
# 
#    Potential and actual number of data values
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "CERI651ALS"

    if action == "setup"
        pbm          = PBM(name)
        pb           = PB(name)
        pb.sifpbname = "CERI651ALS"
        nargin       = length(args)
        pbm.call     = eval( Meta.parse( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["MPOT"] = 10186
        v_["M"] = 61
        v_["MLOWER"] = 9077
        v_["MUPPER"] = 9137
        v_["N"] = 7
        v_["1"] = 1
        v_["ONE"] = 1.0
        v_["X9077"] = 36850.62500
        v_["X9078"] = 36858.00000
        v_["X9079"] = 36865.37500
        v_["X9080"] = 36872.75000
        v_["X9081"] = 36880.12500
        v_["X9082"] = 36887.50000
        v_["X9083"] = 36894.87500
        v_["X9084"] = 36902.25000
        v_["X9085"] = 36909.62500
        v_["X9086"] = 36917.00000
        v_["X9087"] = 36924.37500
        v_["X9088"] = 36931.75000
        v_["X9089"] = 36939.12500
        v_["X9090"] = 36946.50000
        v_["X9091"] = 36953.87500
        v_["X9092"] = 36961.26563
        v_["X9093"] = 36968.67188
        v_["X9094"] = 36976.07813
        v_["X9095"] = 36983.48438
        v_["X9096"] = 36990.89063
        v_["X9097"] = 36998.29688
        v_["X9098"] = 37005.70313
        v_["X9099"] = 37013.10938
        v_["X9100"] = 37020.51563
        v_["X9101"] = 37027.92188
        v_["X9102"] = 37035.32813
        v_["X9103"] = 37042.73438
        v_["X9104"] = 37050.14063
        v_["X9105"] = 37057.54688
        v_["X9106"] = 37064.95313
        v_["X9107"] = 37072.35938
        v_["X9108"] = 37079.76563
        v_["X9109"] = 37087.17188
        v_["X9110"] = 37094.57813
        v_["X9111"] = 37101.98438
        v_["X9112"] = 37109.39063
        v_["X9113"] = 37116.81250
        v_["X9114"] = 37124.25000
        v_["X9115"] = 37131.68750
        v_["X9116"] = 37139.12500
        v_["X9117"] = 37146.56250
        v_["X9118"] = 37154.00000
        v_["X9119"] = 37161.43750
        v_["X9120"] = 37168.87500
        v_["X9121"] = 37176.31250
        v_["X9122"] = 37183.75000
        v_["X9123"] = 37191.18750
        v_["X9124"] = 37198.62500
        v_["X9125"] = 37206.06250
        v_["X9126"] = 37213.50000
        v_["X9127"] = 37220.93750
        v_["X9128"] = 37228.37500
        v_["X9129"] = 37235.81250
        v_["X9130"] = 37243.25000
        v_["X9131"] = 37250.68750
        v_["X9132"] = 37258.12500
        v_["X9133"] = 37265.56250
        v_["X9134"] = 37273.01563
        v_["X9135"] = 37280.48438
        v_["X9136"] = 37287.95313
        v_["X9137"] = 37295.42188
        v_["Y9077"] = 0.00000000
        v_["Y9078"] = 1.96083316
        v_["Y9079"] = 2.94124974
        v_["Y9080"] = 0.98041658
        v_["Y9081"] = 5.88249948
        v_["Y9082"] = 1.96083316
        v_["Y9083"] = 3.92166632
        v_["Y9084"] = 3.92166632
        v_["Y9085"] = 3.92166632
        v_["Y9086"] = 4.90208290
        v_["Y9087"] = 2.94124974
        v_["Y9088"] = 14.70624870
        v_["Y9089"] = 15.68666528
        v_["Y9090"] = 21.56916476
        v_["Y9091"] = 41.17749637
        v_["Y9092"] = 64.70749429
        v_["Y9093"] = 108.82624040
        v_["Y9094"] = 132.35623832
        v_["Y9095"] = 173.53373469
        v_["Y9096"] = 186.27915023
        v_["Y9097"] = 224.51539686
        v_["Y9098"] = 269.61455955
        v_["Y9099"] = 256.86914400
        v_["Y9100"] = 268.63414297
        v_["Y9101"] = 293.14455747
        v_["Y9102"] = 277.45789219
        v_["Y9103"] = 211.76998132
        v_["Y9104"] = 210.78956474
        v_["Y9105"] = 176.47498443
        v_["Y9106"] = 151.96456993
        v_["Y9107"] = 126.47373884
        v_["Y9108"] = 80.39415957
        v_["Y9109"] = 95.10040828
        v_["Y9110"] = 71.57041035
        v_["Y9111"] = 65.68791087
        v_["Y9112"] = 37.25583005
        v_["Y9113"] = 40.19707979
        v_["Y9114"] = 25.49083108
        v_["Y9115"] = 22.54958134
        v_["Y9116"] = 26.47124766
        v_["Y9117"] = 19.60833160
        v_["Y9118"] = 20.58874818
        v_["Y9119"] = 14.70624870
        v_["Y9120"] = 11.76499896
        v_["Y9121"] = 6.86291606
        v_["Y9122"] = 4.90208290
        v_["Y9123"] = 1.96083316
        v_["Y9124"] = 6.86291606
        v_["Y9125"] = 8.82374922
        v_["Y9126"] = 0.98041658
        v_["Y9127"] = 1.96083316
        v_["Y9128"] = 3.92166632
        v_["Y9129"] = 5.88249948
        v_["Y9130"] = 7.84333264
        v_["Y9131"] = 3.92166632
        v_["Y9132"] = 3.92166632
        v_["Y9133"] = 3.92166632
        v_["Y9134"] = 2.94124974
        v_["Y9135"] = 0.98041658
        v_["Y9136"] = 0.98041658
        v_["Y9137"] = 2.94124974
        v_["E9077"] = 1.00000000
        v_["E9078"] = 1.41421356
        v_["E9079"] = 1.73205081
        v_["E9080"] = 1.00000000
        v_["E9081"] = 2.44948974
        v_["E9082"] = 1.41421356
        v_["E9083"] = 2.00000000
        v_["E9084"] = 2.00000000
        v_["E9085"] = 2.00000000
        v_["E9086"] = 2.23606798
        v_["E9087"] = 1.73205081
        v_["E9088"] = 3.87298335
        v_["E9089"] = 4.00000000
        v_["E9090"] = 4.69041576
        v_["E9091"] = 6.48074070
        v_["E9092"] = 8.12403840
        v_["E9093"] = 0.53565375
        v_["E9094"] = 1.61895004
        v_["E9095"] = 3.30413470
        v_["E9096"] = 3.78404875
        v_["E9097"] = 5.13274595
        v_["E9098"] = 6.58312395
        v_["E9099"] = 6.18641406
        v_["E9100"] = 6.55294536
        v_["E9101"] = 7.29161647
        v_["E9102"] = 6.82260384
        v_["E9103"] = 4.69693846
        v_["E9104"] = 4.66287830
        v_["E9105"] = 3.41640786
        v_["E9106"] = 2.44989960
        v_["E9107"] = 1.35781669
        v_["E9108"] = 9.05538514
        v_["E9109"] = 9.84885780
        v_["E9110"] = 8.54400375
        v_["E9111"] = 8.18535277
        v_["E9112"] = 6.16441400
        v_["E9113"] = 6.40312424
        v_["E9114"] = 5.09901951
        v_["E9115"] = 4.79583152
        v_["E9116"] = 5.19615242
        v_["E9117"] = 4.47213595
        v_["E9118"] = 4.58257569
        v_["E9119"] = 3.87298335
        v_["E9120"] = 3.46410162
        v_["E9121"] = 2.64575131
        v_["E9122"] = 2.23606798
        v_["E9123"] = 1.41421356
        v_["E9124"] = 2.64575131
        v_["E9125"] = 3.00000000
        v_["E9126"] = 1.00000000
        v_["E9127"] = 1.41421356
        v_["E9128"] = 2.00000000
        v_["E9129"] = 2.44948974
        v_["E9130"] = 2.82842712
        v_["E9131"] = 2.00000000
        v_["E9132"] = 2.00000000
        v_["E9133"] = 2.00000000
        v_["E9134"] = 1.73205081
        v_["E9135"] = 1.00000000
        v_["E9136"] = 1.00000000
        v_["E9137"] = 1.73205081
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        xscale  = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        iv,ix_,_ = s2x_ii("C",ix_)
        arrset(pb.xnames,iv,"C")
        iv,ix_,_ = s2x_ii("L",ix_)
        arrset(pb.xnames,iv,"L")
        iv,ix_,_ = s2x_ii("A",ix_)
        arrset(pb.xnames,iv,"A")
        iv,ix_,_ = s2x_ii("B",ix_)
        arrset(pb.xnames,iv,"B")
        iv,ix_,_ = s2x_ii("I",ix_)
        arrset(pb.xnames,iv,"I")
        iv,ix_,_ = s2x_ii("S",ix_)
        arrset(pb.xnames,iv,"S")
        iv,ix_,_ = s2x_ii("X0",ix_)
        arrset(pb.xnames,iv,"X0")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["MLOWER"]):Int64(v_["MUPPER"])
            v_["E"] = v_["E"*string(I)]
            v_["EINV"] = v_["ONE"]/v_["E"]
            v_["XOVERE"] = v_["EINV"]*v_["X"*string(I)]
            ig,ig_,_ = s2x_ii("F"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["C"]
            pbm.A[ig,iv] += v_["EINV"]
            iv = ix_["L"]
            pbm.A[ig,iv] += v_["XOVERE"]
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        for I = Int64(v_["MLOWER"]):Int64(v_["MUPPER"])
            v_["E"] = v_["E"*string(I)]
            v_["EINV"] = v_["ONE"]/v_["E"]
            v_["YOVERE"] = v_["EINV"]*v_["Y"*string(I)]
            pbm.gconst[ig_["F"*string(I)]] = v_["YOVERE"]
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.x0[ix_["C"]] = 0.0
        pb.x0[ix_["L"]] = 0.0
        pb.x0[ix_["A"]] = 1.0
        pb.x0[ix_["B"]] = 0.05
        pb.x0[ix_["I"]] = 26061.4
        pb.x0[ix_["S"]] = 38.7105
        pb.x0[ix_["X0"]] = 37027.1
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2x_ii( "B2BEXP", iet_)
        loaset(elftv,it,1,"A")
        loaset(elftv,it,2,"B")
        loaset(elftv,it,3,"I")
        loaset(elftv,it,4,"S")
        loaset(elftv,it,5,"Y")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"X")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["MLOWER"]):Int64(v_["MUPPER"])
            ename = "B"*string(I)
            ie,ie_,_  = s2x_ii(ename,ie_)
            arrset(pbm.elftype,ie,"B2BEXP")
            arrset(ielftype, ie, iet_["B2BEXP"])
            vname = "A"
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="A",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "B"
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="B",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = I
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="I",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "S"
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="S",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X0"
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="X",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,v_["X"*string(I)])
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2x_ii("L2",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for ig in 1:ngrp
            arrset(pbm.grftype,ig,"L2")
        end
        for I = Int64(v_["MLOWER"]):Int64(v_["MUPPER"])
            v_["E"] = v_["E"*string(I)]
            v_["EINV"] = v_["ONE"]/v_["E"]
            ig = ig_["F"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["B"*string(I)])
            loaset(pbm.grelw,ig,posel,v_["EINV"])
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTIONS %%%%%%%
        pb.pbclass = "SUR2-MN-7-0"
        return pb, pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "e_globs"

        pbm = args[1]
        arrset(pbm.efpar,1,sqrt(1.0e0/atan(1.0e0)))    # this is TORPI
        arrset(pbm.efpar,2,sqrt(0.5e0))    # this is ROOTP5
        return pbm

    elseif action == "B2BEXP"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        APB = EV_[1]+EV_[2]
        APB2 = APB*APB
        A2 = EV_[1]*EV_[1]
        B2 = EV_[2]*EV_[2]
        AB = EV_[1]*EV_[2]
        S2 = EV_[4]*EV_[4]
        S3 = EV_[4]*S2
        PI = 0.5e0*AB/APB
        P = EV_[3]*PI
        PAI = 0.5e0*EV_[2]/APB-0.5e0*AB/APB2
        PBI = 0.5e0*EV_[1]/APB-0.5e0*AB/APB2
        PA = PAI*EV_[3]
        loc_PB = PBI*EV_[3]
        PAB = EV_[3]*AB/APB^3
        PAA = -EV_[3]*EV_[2]/APB2+PAB
        PBB = -EV_[3]*EV_[1]/APB2+PAB
        XMY = pbm.elpar[iel_][1]-EV_[5]
        Z = XMY/EV_[4]
        ZY = -1.0e0/EV_[4]
        ZS = -XMY/EV_[4]^2
        ZSY = 1.0e0/EV_[4]^2
        ZSS = 2.0e0*XMY/EV_[4]^3
        R = exp(-0.5e0*Z^2)
        DR = -Z*R
        D2R = -R-Z*DR
        RS = DR*ZS
        RY = DR*ZY
        RSS = D2R*ZS*ZS+DR*ZSS
        RSY = D2R*ZS*ZY+DR*ZSY
        RYY = D2R*ZY*ZY
        AC = pbm.efpar[2]*(EV_[1]*EV_[4]+XMY/EV_[4])
        ACA = pbm.efpar[2]*EV_[4]
        ACS = pbm.efpar[2]*(EV_[1]-XMY/S2)
        ACY = -pbm.efpar[2]/EV_[4]
        ACAS = pbm.efpar[2]
        ACSS = 2.0e0*pbm.efpar[2]*XMY/S3
        ACSY = pbm.efpar[2]/S2
        BC = pbm.efpar[2]*(EV_[2]*EV_[4]+XMY/EV_[4])
        BCB = ACA
        BCS = pbm.efpar[2]*(EV_[2]-XMY/S2)
        BCY = ACY
        BCBS = pbm.efpar[2]
        BCSS = ACSS
        BCSY = ACSY
        QA = ERFC_EV_[4]CEV_[1]LED(AC)
        DQA = 2.0e0*AC*QA-pbm.efpar[1]
        D2QA = 2.0e0*(QA+AC*DQA)
        QAA = DQA*ACA
        QAS = DQA*ACS
        QAY = DQA*ACY
        QAAA = D2QA*ACA*ACA
        QAAS = D2QA*ACA*ACS+DQA*ACAS
        QAAY = D2QA*ACA*ACY
        QASS = D2QA*ACS*ACS+DQA*ACSS
        QASY = D2QA*ACS*ACY+DQA*ACSY
        QAYY = D2QA*ACY*ACY
        QB = ERFC_EV_[4]CEV_[1]LED(BC)
        DQB = 2.0e0*BC*QB-pbm.efpar[1]
        D2QB = 2.0e0*(QB+BC*DQB)
        QBB = DQB*BCB
        QBS = DQB*BCS
        QBY = DQB*BCY
        QBBB = D2QB*BCB*BCB
        QBBS = D2QB*BCB*BCS+DQB*BCBS
        QBBY = D2QB*BCB*BCY
        QBSS = D2QB*BCS*BCS+DQB*BCSS
        QBSY = D2QB*BCS*BCY+DQB*BCSY
        QBYY = D2QB*BCY*BCY
        T = QA+QB
        TA = QAA
        TB = QBB
        TS = QAS+QBS
        TY = QAY+QBY
        TAA = QAAA
        TAS = QAAS
        TAY = QAAY
        TBB = QBBB
        TBS = QBBS
        TBY = QBBY
        TSS = QASS+QBSS
        TSY = QASY+QBSY
        TYY = QAYY+QBYY
        f_   = P*T*R
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = (P*TA+PA*T)*R
            g_[2] = (P*TB+loc_PB*T)*R
            g_[3] = PI*T*R
            g_[4] = P*(T*RS+TS*R)
            g_[5] = P*(T*RY+TY*R)
            if nargout>2
                H_ = zeros(Float64,5,5)
                H_[1,1] = R*(P*TAA+PAA*T+2.0e0*PA*TA)
                H_[1,2] = R*(PA*TB+loc_PB*TA+PAB*T)
                H_[2,1] = H_[1,2]
                H_[1,3] = (PI*TA+PAI*T)*R
                H_[3,1] = H_[1,3]
                H_[1,4] = (P*TA+PA*T)*RS+(P*TAS+PA*TS)*R
                H_[4,1] = H_[1,4]
                H_[1,5] = (P*TA+PA*T)*RY+(P*TAY+PA*TY)*R
                H_[5,1] = H_[1,5]
                H_[2,2] = R*(P*TBB+PBB*T+2.0e0*loc_PB*TB)
                H_[2,3] = (PI*TB+PBI*T)*R
                H_[3,2] = H_[2,3]
                H_[2,4] = (P*TB+loc_PB*T)*RS+(P*TBS+loc_PB*TS)*R
                H_[4,2] = H_[2,4]
                H_[2,5] = (P*TB+loc_PB*T)*RY+(P*TBY+loc_PB*TY)*R
                H_[5,2] = H_[2,5]
                H_[3,4] = PI*(T*RS+TS*R)
                H_[4,3] = H_[3,4]
                H_[3,5] = PI*(T*RY+TY*R)
                H_[5,3] = H_[3,5]
                H_[4,4] = P*(T*RSS+TSS*R+2.0e0*TS*RS)
                H_[4,5] = P*(T*RSY+TSY*R+TS*RY+TY*RS)
                H_[5,4] = H_[4,5]
                H_[5,5] = P*(T*RYY+TYY*R+2.0e0*TY*RY)
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

    elseif action == "L2"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_= GVAR_*GVAR_
        if nargout>1
            g_ = GVAR_+GVAR_
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = 2.0e0
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

    elseif action in  ["fx","fgx","fgHx","cx","cJx","cJHx","cIx","cIJx","cIJHx","cIJxv","fHxv","cJxv","Lxy","Lgxy","LgHxy","LIxy","LIgxy","LIgHxy","LHxyv","LIHxyv"]

        pbm = args[1]
        if pbm.name == name
            pbm.has_globs = [2,0]
            return s2x_eval(action,args...)
        else
            println("ERROR: please run "*name*" with action = setup")
            return ntuple(i->undef,args[end])
        end

    else
        println("ERROR: unknown action "*action*" requested from "*name*"%s.jl")
        return ntuple(i->undef,args[end])
    end

end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

