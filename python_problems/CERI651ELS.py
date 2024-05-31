from s2xlib import *
class  CERI651ELS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

    name = 'CERI651ELS'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'CERI651ELS'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['MPOT'] = 10186
        v_['M'] = 64
        v_['MLOWER'] = 4077
        v_['MUPPER'] = 4140
        v_['N'] = 7
        v_['1'] = 1
        v_['ONE'] = 1.0
        v_['X4077'] = 13558.04688
        v_['X4078'] = 13560.76563
        v_['X4079'] = 13563.48438
        v_['X4080'] = 13566.20313
        v_['X4081'] = 13568.92188
        v_['X4082'] = 13571.64063
        v_['X4083'] = 13574.35938
        v_['X4084'] = 13577.07813
        v_['X4085'] = 13579.79688
        v_['X4086'] = 13582.51563
        v_['X4087'] = 13585.23438
        v_['X4088'] = 13587.95313
        v_['X4089'] = 13590.67188
        v_['X4090'] = 13593.39063
        v_['X4091'] = 13596.10938
        v_['X4092'] = 13598.82813
        v_['X4093'] = 13601.54688
        v_['X4094'] = 13604.26563
        v_['X4095'] = 13606.98438
        v_['X4096'] = 13609.70313
        v_['X4097'] = 13612.42188
        v_['X4098'] = 13615.14063
        v_['X4099'] = 13617.85938
        v_['X4100'] = 13620.57813
        v_['X4101'] = 13623.29688
        v_['X4102'] = 13626.01563
        v_['X4103'] = 13628.73438
        v_['X4104'] = 13631.45313
        v_['X4105'] = 13634.17188
        v_['X4106'] = 13636.89063
        v_['X4107'] = 13639.60938
        v_['X4108'] = 13642.32813
        v_['X4109'] = 13645.04688
        v_['X4110'] = 13647.76563
        v_['X4111'] = 13650.48438
        v_['X4112'] = 13653.20313
        v_['X4113'] = 13655.92188
        v_['X4114'] = 13658.64063
        v_['X4115'] = 13661.35938
        v_['X4116'] = 13664.07813
        v_['X4117'] = 13666.79688
        v_['X4118'] = 13669.51563
        v_['X4119'] = 13672.23438
        v_['X4120'] = 13674.96875
        v_['X4121'] = 13677.71875
        v_['X4122'] = 13680.46875
        v_['X4123'] = 13683.21875
        v_['X4124'] = 13685.96875
        v_['X4125'] = 13688.71875
        v_['X4126'] = 13691.46875
        v_['X4127'] = 13694.21875
        v_['X4128'] = 13696.96875
        v_['X4129'] = 13699.71875
        v_['X4130'] = 13702.46875
        v_['X4131'] = 13705.21875
        v_['X4132'] = 13707.96875
        v_['X4133'] = 13710.71875
        v_['X4134'] = 13713.46875
        v_['X4135'] = 13716.21875
        v_['X4136'] = 13718.96875
        v_['X4137'] = 13721.71875
        v_['X4138'] = 13724.46875
        v_['X4139'] = 13727.21875
        v_['X4140'] = 13729.96875
        v_['Y4077'] = 0.00000000
        v_['Y4078'] = 1.96083316
        v_['Y4079'] = 0.98041658
        v_['Y4080'] = 0.00000000
        v_['Y4081'] = 0.00000000
        v_['Y4082'] = 0.00000000
        v_['Y4083'] = 0.00000000
        v_['Y4084'] = 0.00000000
        v_['Y4085'] = 0.00000000
        v_['Y4086'] = 0.00000000
        v_['Y4087'] = 0.00000000
        v_['Y4088'] = 0.00000000
        v_['Y4089'] = 0.00000000
        v_['Y4090'] = 0.00000000
        v_['Y4091'] = 0.00000000
        v_['Y4092'] = 0.00000000
        v_['Y4093'] = 0.00000000
        v_['Y4094'] = 0.00000000
        v_['Y4095'] = 0.98041658
        v_['Y4096'] = 0.00000000
        v_['Y4097'] = 0.00000000
        v_['Y4098'] = 0.98041658
        v_['Y4099'] = 0.98041658
        v_['Y4100'] = 1.96083316
        v_['Y4101'] = 1.96083316
        v_['Y4102'] = 4.90208290
        v_['Y4103'] = 0.98041658
        v_['Y4104'] = 1.96083316
        v_['Y4105'] = 0.00000000
        v_['Y4106'] = 1.96083316
        v_['Y4107'] = 0.98041658
        v_['Y4108'] = 5.88249948
        v_['Y4109'] = 0.98041658
        v_['Y4110'] = 1.96083316
        v_['Y4111'] = 0.00000000
        v_['Y4112'] = 0.98041658
        v_['Y4113'] = 0.00000000
        v_['Y4114'] = 0.00000000
        v_['Y4115'] = 0.98041658
        v_['Y4116'] = 0.00000000
        v_['Y4117'] = 1.96083316
        v_['Y4118'] = 0.98041658
        v_['Y4119'] = 0.00000000
        v_['Y4120'] = 0.98041658
        v_['Y4121'] = 0.98041658
        v_['Y4122'] = 0.98041658
        v_['Y4123'] = 0.00000000
        v_['Y4124'] = 0.00000000
        v_['Y4125'] = 0.98041658
        v_['Y4126'] = 0.00000000
        v_['Y4127'] = 0.00000000
        v_['Y4128'] = 0.98041658
        v_['Y4129'] = 0.00000000
        v_['Y4130'] = 0.00000000
        v_['Y4131'] = 0.00000000
        v_['Y4132'] = 0.98041658
        v_['Y4133'] = 0.98041658
        v_['Y4134'] = 0.98041658
        v_['Y4135'] = 0.00000000
        v_['Y4136'] = 0.98041658
        v_['Y4137'] = 0.00000000
        v_['Y4138'] = 1.96083316
        v_['Y4139'] = 0.00000000
        v_['Y4140'] = 0.00000000
        v_['E4077'] = 1.00000000
        v_['E4078'] = 1.41421356
        v_['E4079'] = 1.00000000
        v_['E4080'] = 1.00000000
        v_['E4081'] = 1.00000000
        v_['E4082'] = 1.00000000
        v_['E4083'] = 1.00000000
        v_['E4084'] = 1.00000000
        v_['E4085'] = 1.00000000
        v_['E4086'] = 1.00000000
        v_['E4087'] = 1.00000000
        v_['E4088'] = 1.00000000
        v_['E4089'] = 1.00000000
        v_['E4090'] = 1.00000000
        v_['E4091'] = 1.00000000
        v_['E4092'] = 1.00000000
        v_['E4093'] = 1.00000000
        v_['E4094'] = 1.00000000
        v_['E4095'] = 1.00000000
        v_['E4096'] = 1.00000000
        v_['E4097'] = 1.00000000
        v_['E4098'] = 1.00000000
        v_['E4099'] = 1.00000000
        v_['E4100'] = 1.41421356
        v_['E4101'] = 1.41421356
        v_['E4102'] = 2.23606798
        v_['E4103'] = 1.00000000
        v_['E4104'] = 1.41421356
        v_['E4105'] = 1.00000000
        v_['E4106'] = 1.41421356
        v_['E4107'] = 1.00000000
        v_['E4108'] = 2.44948974
        v_['E4109'] = 1.00000000
        v_['E4110'] = 1.41421356
        v_['E4111'] = 1.00000000
        v_['E4112'] = 1.00000000
        v_['E4113'] = 1.00000000
        v_['E4114'] = 1.00000000
        v_['E4115'] = 1.00000000
        v_['E4116'] = 1.00000000
        v_['E4117'] = 1.41421356
        v_['E4118'] = 1.00000000
        v_['E4119'] = 1.00000000
        v_['E4120'] = 1.00000000
        v_['E4121'] = 1.00000000
        v_['E4122'] = 1.00000000
        v_['E4123'] = 1.00000000
        v_['E4124'] = 1.00000000
        v_['E4125'] = 1.00000000
        v_['E4126'] = 1.00000000
        v_['E4127'] = 1.00000000
        v_['E4128'] = 1.00000000
        v_['E4129'] = 1.00000000
        v_['E4130'] = 1.00000000
        v_['E4131'] = 1.00000000
        v_['E4132'] = 1.00000000
        v_['E4133'] = 1.00000000
        v_['E4134'] = 1.00000000
        v_['E4135'] = 1.00000000
        v_['E4136'] = 1.00000000
        v_['E4137'] = 1.00000000
        v_['E4138'] = 1.41421356
        v_['E4139'] = 1.00000000
        v_['E4140'] = 1.00000000
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
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['C']
            pbm.A[ig,iv] = v_['EINV']+pbm.A[ig,iv]
            iv = ix_['L']
            pbm.A[ig,iv] = v_['XOVERE']+pbm.A[ig,iv]
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
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
        pb.x0[ix_['C']] = 0.0
        pb.x0[ix_['L']] = 0.0
        pb.x0[ix_['A']] = 1.0
        pb.x0[ix_['B']] = 0.05
        pb.x0[ix_['I']] = 17.06794
        pb.x0[ix_['S']] = 8.0
        pb.x0[ix_['X0']] = 13642.3
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
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2x_ii('L2',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for ig in range(0,ngrp):
            pbm.grftype = arrset(pbm.grftype,ig,'L2')
        for I in range(int(v_['MLOWER']),int(v_['MUPPER'])+1):
            v_['E'] = v_['E'+str(I)]
            v_['EINV'] = v_['ONE']/v_['E']
            ig = ig_['F'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['B'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel,v_['EINV'])
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "SUR2-MN-7-0"
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

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def L2(pbm,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= GVAR_*GVAR_
        if nargout>1:
            g_ = GVAR_+GVAR_
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 2.0e0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

