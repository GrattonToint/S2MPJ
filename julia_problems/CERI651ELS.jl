function CERI651ELS(action,args...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CERI651ELS
#    *********
# 
#    ISIS Data fitting problem CERI651E given as an inconsistent set of
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
#    subset X in [13556.2988352, 13731.2988352]
# 
#    SIF input: Nick Gould and Tyrone Rees, Mar 2016
#    Least-squares version of CERI651E.SIF, Nick Gould, Jan 2020.
# 
#    classification = "SUR2-MN-7-0"
# 
#    Potential and actual number of data values
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "CERI651ELS"

    if action == "setup"
        pbm          = PBM(name)
        pb           = PB(name)
        pb.sifpbname = "CERI651ELS"
        nargin       = length(args)
        pbm.call     = eval( Meta.parse( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["MPOT"] = 10186
        v_["M"] = 64
        v_["MLOWER"] = 4077
        v_["MUPPER"] = 4140
        v_["N"] = 7
        v_["1"] = 1
        v_["ONE"] = 1.0
        v_["X4077"] = 13558.04688
        v_["X4078"] = 13560.76563
        v_["X4079"] = 13563.48438
        v_["X4080"] = 13566.20313
        v_["X4081"] = 13568.92188
        v_["X4082"] = 13571.64063
        v_["X4083"] = 13574.35938
        v_["X4084"] = 13577.07813
        v_["X4085"] = 13579.79688
        v_["X4086"] = 13582.51563
        v_["X4087"] = 13585.23438
        v_["X4088"] = 13587.95313
        v_["X4089"] = 13590.67188
        v_["X4090"] = 13593.39063
        v_["X4091"] = 13596.10938
        v_["X4092"] = 13598.82813
        v_["X4093"] = 13601.54688
        v_["X4094"] = 13604.26563
        v_["X4095"] = 13606.98438
        v_["X4096"] = 13609.70313
        v_["X4097"] = 13612.42188
        v_["X4098"] = 13615.14063
        v_["X4099"] = 13617.85938
        v_["X4100"] = 13620.57813
        v_["X4101"] = 13623.29688
        v_["X4102"] = 13626.01563
        v_["X4103"] = 13628.73438
        v_["X4104"] = 13631.45313
        v_["X4105"] = 13634.17188
        v_["X4106"] = 13636.89063
        v_["X4107"] = 13639.60938
        v_["X4108"] = 13642.32813
        v_["X4109"] = 13645.04688
        v_["X4110"] = 13647.76563
        v_["X4111"] = 13650.48438
        v_["X4112"] = 13653.20313
        v_["X4113"] = 13655.92188
        v_["X4114"] = 13658.64063
        v_["X4115"] = 13661.35938
        v_["X4116"] = 13664.07813
        v_["X4117"] = 13666.79688
        v_["X4118"] = 13669.51563
        v_["X4119"] = 13672.23438
        v_["X4120"] = 13674.96875
        v_["X4121"] = 13677.71875
        v_["X4122"] = 13680.46875
        v_["X4123"] = 13683.21875
        v_["X4124"] = 13685.96875
        v_["X4125"] = 13688.71875
        v_["X4126"] = 13691.46875
        v_["X4127"] = 13694.21875
        v_["X4128"] = 13696.96875
        v_["X4129"] = 13699.71875
        v_["X4130"] = 13702.46875
        v_["X4131"] = 13705.21875
        v_["X4132"] = 13707.96875
        v_["X4133"] = 13710.71875
        v_["X4134"] = 13713.46875
        v_["X4135"] = 13716.21875
        v_["X4136"] = 13718.96875
        v_["X4137"] = 13721.71875
        v_["X4138"] = 13724.46875
        v_["X4139"] = 13727.21875
        v_["X4140"] = 13729.96875
        v_["Y4077"] = 0.00000000
        v_["Y4078"] = 1.96083316
        v_["Y4079"] = 0.98041658
        v_["Y4080"] = 0.00000000
        v_["Y4081"] = 0.00000000
        v_["Y4082"] = 0.00000000
        v_["Y4083"] = 0.00000000
        v_["Y4084"] = 0.00000000
        v_["Y4085"] = 0.00000000
        v_["Y4086"] = 0.00000000
        v_["Y4087"] = 0.00000000
        v_["Y4088"] = 0.00000000
        v_["Y4089"] = 0.00000000
        v_["Y4090"] = 0.00000000
        v_["Y4091"] = 0.00000000
        v_["Y4092"] = 0.00000000
        v_["Y4093"] = 0.00000000
        v_["Y4094"] = 0.00000000
        v_["Y4095"] = 0.98041658
        v_["Y4096"] = 0.00000000
        v_["Y4097"] = 0.00000000
        v_["Y4098"] = 0.98041658
        v_["Y4099"] = 0.98041658
        v_["Y4100"] = 1.96083316
        v_["Y4101"] = 1.96083316
        v_["Y4102"] = 4.90208290
        v_["Y4103"] = 0.98041658
        v_["Y4104"] = 1.96083316
        v_["Y4105"] = 0.00000000
        v_["Y4106"] = 1.96083316
        v_["Y4107"] = 0.98041658
        v_["Y4108"] = 5.88249948
        v_["Y4109"] = 0.98041658
        v_["Y4110"] = 1.96083316
        v_["Y4111"] = 0.00000000
        v_["Y4112"] = 0.98041658
        v_["Y4113"] = 0.00000000
        v_["Y4114"] = 0.00000000
        v_["Y4115"] = 0.98041658
        v_["Y4116"] = 0.00000000
        v_["Y4117"] = 1.96083316
        v_["Y4118"] = 0.98041658
        v_["Y4119"] = 0.00000000
        v_["Y4120"] = 0.98041658
        v_["Y4121"] = 0.98041658
        v_["Y4122"] = 0.98041658
        v_["Y4123"] = 0.00000000
        v_["Y4124"] = 0.00000000
        v_["Y4125"] = 0.98041658
        v_["Y4126"] = 0.00000000
        v_["Y4127"] = 0.00000000
        v_["Y4128"] = 0.98041658
        v_["Y4129"] = 0.00000000
        v_["Y4130"] = 0.00000000
        v_["Y4131"] = 0.00000000
        v_["Y4132"] = 0.98041658
        v_["Y4133"] = 0.98041658
        v_["Y4134"] = 0.98041658
        v_["Y4135"] = 0.00000000
        v_["Y4136"] = 0.98041658
        v_["Y4137"] = 0.00000000
        v_["Y4138"] = 1.96083316
        v_["Y4139"] = 0.00000000
        v_["Y4140"] = 0.00000000
        v_["E4077"] = 1.00000000
        v_["E4078"] = 1.41421356
        v_["E4079"] = 1.00000000
        v_["E4080"] = 1.00000000
        v_["E4081"] = 1.00000000
        v_["E4082"] = 1.00000000
        v_["E4083"] = 1.00000000
        v_["E4084"] = 1.00000000
        v_["E4085"] = 1.00000000
        v_["E4086"] = 1.00000000
        v_["E4087"] = 1.00000000
        v_["E4088"] = 1.00000000
        v_["E4089"] = 1.00000000
        v_["E4090"] = 1.00000000
        v_["E4091"] = 1.00000000
        v_["E4092"] = 1.00000000
        v_["E4093"] = 1.00000000
        v_["E4094"] = 1.00000000
        v_["E4095"] = 1.00000000
        v_["E4096"] = 1.00000000
        v_["E4097"] = 1.00000000
        v_["E4098"] = 1.00000000
        v_["E4099"] = 1.00000000
        v_["E4100"] = 1.41421356
        v_["E4101"] = 1.41421356
        v_["E4102"] = 2.23606798
        v_["E4103"] = 1.00000000
        v_["E4104"] = 1.41421356
        v_["E4105"] = 1.00000000
        v_["E4106"] = 1.41421356
        v_["E4107"] = 1.00000000
        v_["E4108"] = 2.44948974
        v_["E4109"] = 1.00000000
        v_["E4110"] = 1.41421356
        v_["E4111"] = 1.00000000
        v_["E4112"] = 1.00000000
        v_["E4113"] = 1.00000000
        v_["E4114"] = 1.00000000
        v_["E4115"] = 1.00000000
        v_["E4116"] = 1.00000000
        v_["E4117"] = 1.41421356
        v_["E4118"] = 1.00000000
        v_["E4119"] = 1.00000000
        v_["E4120"] = 1.00000000
        v_["E4121"] = 1.00000000
        v_["E4122"] = 1.00000000
        v_["E4123"] = 1.00000000
        v_["E4124"] = 1.00000000
        v_["E4125"] = 1.00000000
        v_["E4126"] = 1.00000000
        v_["E4127"] = 1.00000000
        v_["E4128"] = 1.00000000
        v_["E4129"] = 1.00000000
        v_["E4130"] = 1.00000000
        v_["E4131"] = 1.00000000
        v_["E4132"] = 1.00000000
        v_["E4133"] = 1.00000000
        v_["E4134"] = 1.00000000
        v_["E4135"] = 1.00000000
        v_["E4136"] = 1.00000000
        v_["E4137"] = 1.00000000
        v_["E4138"] = 1.41421356
        v_["E4139"] = 1.00000000
        v_["E4140"] = 1.00000000
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
        pb.x0[ix_["I"]] = 17.06794
        pb.x0[ix_["S"]] = 8.0
        pb.x0[ix_["X0"]] = 13642.3
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

