from s2xlib import *
class  CERI651A(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CERI651A
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
# 
#    classification = "NOR2-MN-7-61"
# 
#    Potential and actual number of data values
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'CERI651A'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'CERI651A'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['MPOT'] = 10186
        v_['M'] = 61
        v_['MLOWER'] = 9077
        v_['MUPPER'] = 9137
        v_['N'] = 7
        v_['1'] = 1
        v_['ONE'] = 1.0
        v_['X9077'] = 36850.62500
        v_['X9078'] = 36858.00000
        v_['X9079'] = 36865.37500
        v_['X9080'] = 36872.75000
        v_['X9081'] = 36880.12500
        v_['X9082'] = 36887.50000
        v_['X9083'] = 36894.87500
        v_['X9084'] = 36902.25000
        v_['X9085'] = 36909.62500
        v_['X9086'] = 36917.00000
        v_['X9087'] = 36924.37500
        v_['X9088'] = 36931.75000
        v_['X9089'] = 36939.12500
        v_['X9090'] = 36946.50000
        v_['X9091'] = 36953.87500
        v_['X9092'] = 36961.26563
        v_['X9093'] = 36968.67188
        v_['X9094'] = 36976.07813
        v_['X9095'] = 36983.48438
        v_['X9096'] = 36990.89063
        v_['X9097'] = 36998.29688
        v_['X9098'] = 37005.70313
        v_['X9099'] = 37013.10938
        v_['X9100'] = 37020.51563
        v_['X9101'] = 37027.92188
        v_['X9102'] = 37035.32813
        v_['X9103'] = 37042.73438
        v_['X9104'] = 37050.14063
        v_['X9105'] = 37057.54688
        v_['X9106'] = 37064.95313
        v_['X9107'] = 37072.35938
        v_['X9108'] = 37079.76563
        v_['X9109'] = 37087.17188
        v_['X9110'] = 37094.57813
        v_['X9111'] = 37101.98438
        v_['X9112'] = 37109.39063
        v_['X9113'] = 37116.81250
        v_['X9114'] = 37124.25000
        v_['X9115'] = 37131.68750
        v_['X9116'] = 37139.12500
        v_['X9117'] = 37146.56250
        v_['X9118'] = 37154.00000
        v_['X9119'] = 37161.43750
        v_['X9120'] = 37168.87500
        v_['X9121'] = 37176.31250
        v_['X9122'] = 37183.75000
        v_['X9123'] = 37191.18750
        v_['X9124'] = 37198.62500
        v_['X9125'] = 37206.06250
        v_['X9126'] = 37213.50000
        v_['X9127'] = 37220.93750
        v_['X9128'] = 37228.37500
        v_['X9129'] = 37235.81250
        v_['X9130'] = 37243.25000
        v_['X9131'] = 37250.68750
        v_['X9132'] = 37258.12500
        v_['X9133'] = 37265.56250
        v_['X9134'] = 37273.01563
        v_['X9135'] = 37280.48438
        v_['X9136'] = 37287.95313
        v_['X9137'] = 37295.42188
        v_['Y9077'] = 0.00000000
        v_['Y9078'] = 1.96083316
        v_['Y9079'] = 2.94124974
        v_['Y9080'] = 0.98041658
        v_['Y9081'] = 5.88249948
        v_['Y9082'] = 1.96083316
        v_['Y9083'] = 3.92166632
        v_['Y9084'] = 3.92166632
        v_['Y9085'] = 3.92166632
        v_['Y9086'] = 4.90208290
        v_['Y9087'] = 2.94124974
        v_['Y9088'] = 14.70624870
        v_['Y9089'] = 15.68666528
        v_['Y9090'] = 21.56916476
        v_['Y9091'] = 41.17749637
        v_['Y9092'] = 64.70749429
        v_['Y9093'] = 108.82624040
        v_['Y9094'] = 132.35623832
        v_['Y9095'] = 173.53373469
        v_['Y9096'] = 186.27915023
        v_['Y9097'] = 224.51539686
        v_['Y9098'] = 269.61455955
        v_['Y9099'] = 256.86914400
        v_['Y9100'] = 268.63414297
        v_['Y9101'] = 293.14455747
        v_['Y9102'] = 277.45789219
        v_['Y9103'] = 211.76998132
        v_['Y9104'] = 210.78956474
        v_['Y9105'] = 176.47498443
        v_['Y9106'] = 151.96456993
        v_['Y9107'] = 126.47373884
        v_['Y9108'] = 80.39415957
        v_['Y9109'] = 95.10040828
        v_['Y9110'] = 71.57041035
        v_['Y9111'] = 65.68791087
        v_['Y9112'] = 37.25583005
        v_['Y9113'] = 40.19707979
        v_['Y9114'] = 25.49083108
        v_['Y9115'] = 22.54958134
        v_['Y9116'] = 26.47124766
        v_['Y9117'] = 19.60833160
        v_['Y9118'] = 20.58874818
        v_['Y9119'] = 14.70624870
        v_['Y9120'] = 11.76499896
        v_['Y9121'] = 6.86291606
        v_['Y9122'] = 4.90208290
        v_['Y9123'] = 1.96083316
        v_['Y9124'] = 6.86291606
        v_['Y9125'] = 8.82374922
        v_['Y9126'] = 0.98041658
        v_['Y9127'] = 1.96083316
        v_['Y9128'] = 3.92166632
        v_['Y9129'] = 5.88249948
        v_['Y9130'] = 7.84333264
        v_['Y9131'] = 3.92166632
        v_['Y9132'] = 3.92166632
        v_['Y9133'] = 3.92166632
        v_['Y9134'] = 2.94124974
        v_['Y9135'] = 0.98041658
        v_['Y9136'] = 0.98041658
        v_['Y9137'] = 2.94124974
        v_['E9077'] = 1.00000000
        v_['E9078'] = 1.41421356
        v_['E9079'] = 1.73205081
        v_['E9080'] = 1.00000000
        v_['E9081'] = 2.44948974
        v_['E9082'] = 1.41421356
        v_['E9083'] = 2.00000000
        v_['E9084'] = 2.00000000
        v_['E9085'] = 2.00000000
        v_['E9086'] = 2.23606798
        v_['E9087'] = 1.73205081
        v_['E9088'] = 3.87298335
        v_['E9089'] = 4.00000000
        v_['E9090'] = 4.69041576
        v_['E9091'] = 6.48074070
        v_['E9092'] = 8.12403840
        v_['E9093'] = 0.53565375
        v_['E9094'] = 1.61895004
        v_['E9095'] = 3.30413470
        v_['E9096'] = 3.78404875
        v_['E9097'] = 5.13274595
        v_['E9098'] = 6.58312395
        v_['E9099'] = 6.18641406
        v_['E9100'] = 6.55294536
        v_['E9101'] = 7.29161647
        v_['E9102'] = 6.82260384
        v_['E9103'] = 4.69693846
        v_['E9104'] = 4.66287830
        v_['E9105'] = 3.41640786
        v_['E9106'] = 2.44989960
        v_['E9107'] = 1.35781669
        v_['E9108'] = 9.05538514
        v_['E9109'] = 9.84885780
        v_['E9110'] = 8.54400375
        v_['E9111'] = 8.18535277
        v_['E9112'] = 6.16441400
        v_['E9113'] = 6.40312424
        v_['E9114'] = 5.09901951
        v_['E9115'] = 4.79583152
        v_['E9116'] = 5.19615242
        v_['E9117'] = 4.47213595
        v_['E9118'] = 4.58257569
        v_['E9119'] = 3.87298335
        v_['E9120'] = 3.46410162
        v_['E9121'] = 2.64575131
        v_['E9122'] = 2.23606798
        v_['E9123'] = 1.41421356
        v_['E9124'] = 2.64575131
        v_['E9125'] = 3.00000000
        v_['E9126'] = 1.00000000
        v_['E9127'] = 1.41421356
        v_['E9128'] = 2.00000000
        v_['E9129'] = 2.44948974
        v_['E9130'] = 2.82842712
        v_['E9131'] = 2.00000000
        v_['E9132'] = 2.00000000
        v_['E9133'] = 2.00000000
        v_['E9134'] = 1.73205081
        v_['E9135'] = 1.00000000
        v_['E9136'] = 1.00000000
        v_['E9137'] = 1.73205081
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2x_ii('C',ix_)
        pb.xnames=arrset(pb.xnames,iv,'C')
        [iv,ix_,_] = s2x_ii('L',ix_)
        pb.xnames=arrset(pb.xnames,iv,'L')
        [iv,ix_,_] = s2x_ii('A',ix_)
        pb.xnames=arrset(pb.xnames,iv,'A')
        [iv,ix_,_] = s2x_ii('B',ix_)
        pb.xnames=arrset(pb.xnames,iv,'B')
        [iv,ix_,_] = s2x_ii('I',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I')
        [iv,ix_,_] = s2x_ii('S',ix_)
        pb.xnames=arrset(pb.xnames,iv,'S')
        [iv,ix_,_] = s2x_ii('X0',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X0')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['MLOWER']),int(v_['MUPPER'])+1):
            v_['E'] = v_['E'+str(I)]
            v_['EINV'] = v_['ONE']/v_['E']
            v_['XOVERE'] = v_['EINV']*v_['X'+str(I)]
            [ig,ig_,_] = s2x_ii('F'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'F'+str(I))
            iv = ix_['C']
            pbm.A[ig,iv] = v_['EINV']+pbm.A[ig,iv]
            iv = ix_['L']
            pbm.A[ig,iv] = v_['XOVERE']+pbm.A[ig,iv]
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        legrps = find(gtype,lambda x:x=='<=')
        eqgrps = find(gtype,lambda x:x=='==')
        gegrps = find(gtype,lambda x:x=='>=')
        pb.nle = len(legrps)
        pb.neq = len(eqgrps)
        pb.nge = len(gegrps)
        pb.m   = pb.nle+pb.neq+pb.nge
        pbm.congrps = find(gtype,lambda x:(x=='<=' or x=='==' or x=='>='))
        pb.cnames= cnames[pbm.congrps]
        pb.nob = ngrp-pb.m
        pbm.objgrps = find(gtype,lambda x:x=='<>')
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        for I in range(int(v_['MLOWER']),int(v_['MUPPER'])+1):
            v_['E'] = v_['E'+str(I)]
            v_['EINV'] = v_['ONE']/v_['E']
            v_['YOVERE'] = v_['EINV']*v_['Y'+str(I)]
            pbm.gconst = arrset(pbm.gconst,ig_['F'+str(I)],v_['YOVERE'])
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        if('C' in ix_):
            pb.x0[ix_['C']] = 0.0
        else:
            pb.y0 = arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['C']),0.0)
        if('L' in ix_):
            pb.x0[ix_['L']] = 0.0
        else:
            pb.y0 = arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['L']),0.0)
        if('A' in ix_):
            pb.x0[ix_['A']] = 1.0
        else:
            pb.y0 = arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['A']),1.0)
        if('B' in ix_):
            pb.x0[ix_['B']] = 0.05
        else:
            pb.y0 = arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['B']),0.05)
        if('I' in ix_):
            pb.x0[ix_['I']] = 26061.4
        else:
            pb.y0 = arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['I']),26061.4)
        if('S' in ix_):
            pb.x0[ix_['S']] = 38.7105
        else:
            pb.y0 = arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['S']),38.7105)
        if('X0' in ix_):
            pb.x0[ix_['X0']] = 37027.1
        else:
            pb.y0 = arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X0']),37027.1)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'B2BEXP', iet_)
        elftv = loaset(elftv,it,0,'A')
        elftv = loaset(elftv,it,1,'B')
        elftv = loaset(elftv,it,2,'I')
        elftv = loaset(elftv,it,3,'S')
        elftv = loaset(elftv,it,4,'Y')
        elftp = []
        elftp = loaset(elftp,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['MLOWER']),int(v_['MUPPER'])+1):
            ename = 'B'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'B2BEXP')
            ielftype = arrset(ielftype, ie, iet_["B2BEXP"])
            vname = 'A'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='A')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'B'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='B')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = I
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='I')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'S'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='S')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X0'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='X')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],v_['X'+str(I)])
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['MLOWER']),int(v_['MUPPER'])+1):
            v_['E'] = v_['E'+str(I)]
            v_['EINV'] = v_['ONE']/v_['E']
            ig = ig_['F'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['B'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,v_['EINV'])
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "NOR2-MN-7-61"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def e_globs(pbm):

        import numpy as np
        pbm.efpar = np.array([]);
        pbm.efpar  = (
              arrset( pbm.efpar,0,np.sqrt(1.0e0/np.arctan(1.0e0)))    # this is TORPI)
        pbm.efpar = arrset( pbm.efpar,1,np.sqrt(0.5e0))    # this is ROOTP5
        return pbm

    @staticmethod
    def B2BEXP(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        APB = EV_[0]+EV_[1]
        APB2 = APB*APB
        A2 = EV_[0]*EV_[0]
        B2 = EV_[1]*EV_[1]
        AB = EV_[0]*EV_[1]
        S2 = EV_[3]*EV_[3]
        S3 = EV_[3]*S2
        PI = 0.5e0*AB/APB
        P = EV_[2]*PI
        PAI = 0.5e0*EV_[1]/APB-0.5e0*AB/APB2
        PBI = 0.5e0*EV_[0]/APB-0.5e0*AB/APB2
        PA = PAI*EV_[2]
        loc_PB = PBI*EV_[2]
        PAB = EV_[2]*AB/APB**3
        PAA = -EV_[2]*EV_[1]/APB2+PAB
        PBB = -EV_[2]*EV_[0]/APB2+PAB
        XMY = pbm.elpar[iel_][0]-EV_[4]
        Z = XMY/EV_[3]
        ZY = -1.0e0/EV_[3]
        ZS = -XMY/EV_[3]**2
        ZSY = 1.0e0/EV_[3]**2
        ZSS = 2.0e0*XMY/EV_[3]**3
        R = np.exp(-0.5e0*Z**2)
        DR = -Z*R
        D2R = -R-Z*DR
        RS = DR*ZS
        RY = DR*ZY
        RSS = D2R*ZS*ZS+DR*ZSS
        RSY = D2R*ZS*ZY+DR*ZSY
        RYY = D2R*ZY*ZY
        AC = pbm.efpar[1]*(EV_[0]*EV_[3]+XMY/EV_[3])
        ACA = pbm.efpar[1]*EV_[3]
        ACS = pbm.efpar[1]*(EV_[0]-XMY/S2)
        ACY = -pbm.efpar[1]/EV_[3]
        ACAS = pbm.efpar[1]
        ACSS = 2.0e0*pbm.efpar[1]*XMY/S3
        ACSY = pbm.efpar[1]/S2
        BC = pbm.efpar[1]*(EV_[1]*EV_[3]+XMY/EV_[3])
        BCB = ACA
        BCS = pbm.efpar[1]*(EV_[1]-XMY/S2)
        BCY = ACY
        BCBS = pbm.efpar[1]
        BCSS = ACSS
        BCSY = ACSY
        QA = ERFC_EV_[3]CEV_[0]LED(AC)
        DQA = 2.0e0*AC*QA-pbm.efpar[0]
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
        QB = ERFC_EV_[3]CEV_[0]LED(BC)
        DQB = 2.0e0*BC*QB-pbm.efpar[0]
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
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = (P*TA+PA*T)*R
            g_[1] = (P*TB+loc_PB*T)*R
            g_[2] = PI*T*R
            g_[3] = P*(T*RS+TS*R)
            g_[4] = P*(T*RY+TY*R)
            if nargout>2:
                H_ = np.zeros((5,5))
                H_[0,0] = R*(P*TAA+PAA*T+2.0e0*PA*TA)
                H_[0,1] = R*(PA*TB+loc_PB*TA+PAB*T)
                H_[1,0] = H_[0,1]
                H_[0,2] = (PI*TA+PAI*T)*R
                H_[2,0] = H_[0,2]
                H_[0,3] = (P*TA+PA*T)*RS+(P*TAS+PA*TS)*R
                H_[3,0] = H_[0,3]
                H_[0,4] = (P*TA+PA*T)*RY+(P*TAY+PA*TY)*R
                H_[4,0] = H_[0,4]
                H_[1,1] = R*(P*TBB+PBB*T+2.0e0*loc_PB*TB)
                H_[1,2] = (PI*TB+PBI*T)*R
                H_[2,1] = H_[1,2]
                H_[1,3] = (P*TB+loc_PB*T)*RS+(P*TBS+loc_PB*TS)*R
                H_[3,1] = H_[1,3]
                H_[1,4] = (P*TB+loc_PB*T)*RY+(P*TBY+loc_PB*TY)*R
                H_[4,1] = H_[1,4]
                H_[2,3] = PI*(T*RS+TS*R)
                H_[3,2] = H_[2,3]
                H_[2,4] = PI*(T*RY+TY*R)
                H_[4,2] = H_[2,4]
                H_[3,3] = P*(T*RSS+TSS*R+2.0e0*TS*RS)
                H_[3,4] = P*(T*RSY+TSY*R+TS*RY+TY*RS)
                H_[4,3] = H_[3,4]
                H_[4,4] = P*(T*RYY+TYY*R+2.0e0*TY*RY)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

