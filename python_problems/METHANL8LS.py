from s2mpjlib import *
class  METHANL8LS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    The methanol-8 problem by Fletcher.
# 
#    Source: Problem 2c in
#    J.J. More',"A collection of nonlinear model problems"
#    Proceedings of the AMS-SIAM Summer Seminar on the Computational
#    Solution of Nonlinear Systems of Equations, Colorado, 1988.
#    Argonne National Laboratory MCS-P60-0289, 1989.
# 
#    SIF input: N. Gould, Dec 1989.
#    Least-squares version of METHANL8.SIF, Nick Gould, Jan 2020.
# 
#    classification = "C-CSUR2-MN-31-0"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'METHANL8LS'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 8
        v_['M'] = 2
        v_['K'] = 2
        v_['0'] = 0
        v_['1'] = 1
        v_['N-1'] = -1+v_['N']
        v_['N-2'] = -2+v_['N']
        v_['K-'] = -1+v_['K']
        v_['K+'] = 1+v_['K']
        v_['A1'] = 18.5751
        v_['B1'] = -3632.649
        v_['C1'] = 239.2
        v_['A2'] = 18.3443
        v_['B2'] = -3841.2203
        v_['C2'] = 228.0
        v_['AL1'] = 0.0
        v_['ALp1'] = 15.97
        v_['ALpp1'] = 0.0422
        v_['AL2'] = 0.0
        v_['ALp2'] = 18.1
        v_['ALpp2'] = 0.0
        v_['BE1'] = 9566.67
        v_['BEp1'] = -1.59
        v_['BEpp1'] = 0.0422
        v_['BE2'] = 10834.67
        v_['BEp2'] = 8.74
        v_['BEpp2'] = 0.0
        v_['FL1'] = 451.25
        v_['FL2'] = 684.25
        v_['FV1'] = 0.0
        v_['FV2'] = 0.0
        v_['TF'] = 89.0
        v_['B'] = 693.37
        v_['D'] = 442.13
        v_['Q'] = 8386200.0
        v_['PI0'] = 1210.0
        v_['PI1'] = 1200.0
        v_['PI2'] = 1190.0
        v_['PI3'] = 1180.0
        v_['PI4'] = 1170.0
        v_['PI5'] = 1160.0
        v_['PI6'] = 1150.0
        v_['PI7'] = 1140.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['0']),int(v_['N-1'])+1):
            [iv,ix_,_] = s2mpj_ii('T'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'T'+str(I))
            v_['INVPI'+str(I)] = 1.0/v_['PI'+str(I)]
            for J in range(int(v_['1']),int(v_['M'])+1):
                [iv,ix_,_] = s2mpj_ii('X'+str(I)+','+str(J),ix_)
                self.xnames=arrset(self.xnames,iv,'X'+str(I)+','+str(J))
        for I in range(int(v_['0']),int(v_['N-2'])+1):
            [iv,ix_,_] = s2mpj_ii('V'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'V'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for J in range(int(v_['1']),int(v_['M'])+1):
            [ig,ig_,_] = s2mpj_ii('2.1-'+str(J),ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['0']))+','+str(J)]])
            valA = np.append(valA,float(v_['B']))
            self.gscale = arrset(self.gscale,ig,float(1.0e+4))
            [ig,ig_,_] = s2mpj_ii('2.3-'+str(J),ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['N-1']))+','+str(J)]])
            valA = np.append(valA,float(-1.0))
            for I in range(int(v_['1']),int(v_['N-2'])+1):
                [ig,ig_,_] = s2mpj_ii('2.2-'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
                self.gscale = arrset(self.gscale,ig,float(1.0e+4))
        for I in range(int(v_['0']),int(v_['N-1'])+1):
            [ig,ig_,_] = s2mpj_ii('2.7-'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('2.8',ig_)
        gtype = arrset(gtype,ig,'<>')
        self.gscale = arrset(self.gscale,ig,float(1.0e+10))
        for I in range(int(v_['1']),int(v_['N-2'])+1):
            [ig,ig_,_] = s2mpj_ii('2.9-'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            self.gscale = arrset(self.gscale,ig,float(1.0e+10))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        v_['SMALLHF'] = 0.0e+0
        v_['BIGHF'] = 0.0e+0
        for J in range(int(v_['1']),int(v_['M'])+1):
            self.gconst  = (
                  arrset(self.gconst,ig_['2.2-'+str(int(v_['K']))+','+str(J)],float(v_['FL'+str(J)])))
            self.gconst  = (
                  arrset(self.gconst,ig_['2.2-'+str(int(v_['K+']))+','+str(J)],float(v_['FV'+str(J)])))
            v_['TFTF'] = v_['TF']*v_['TF']
            v_['TEMP1'] = v_['TFTF']*v_['ALpp'+str(J)]
            v_['TEMP2'] = v_['TF']*v_['ALp'+str(J)]
            v_['TEMP1'] = v_['TEMP1']+v_['TEMP2']
            v_['TEMP1'] = v_['TEMP1']+v_['AL'+str(J)]
            v_['TEMP1'] = v_['TEMP1']*v_['FL'+str(J)]
            v_['SMALLHF'] = v_['SMALLHF']+v_['TEMP1']
            v_['TEMP1'] = v_['TFTF']*v_['BEpp'+str(J)]
            v_['TEMP2'] = v_['TF']*v_['BEp'+str(J)]
            v_['TEMP1'] = v_['TEMP1']+v_['TEMP2']
            v_['TEMP1'] = v_['TEMP1']+v_['BE'+str(J)]
            v_['TEMP1'] = v_['TEMP1']*v_['FV'+str(J)]
            v_['BIGHF'] = v_['BIGHF']+v_['TEMP1']
        for I in range(int(v_['0']),int(v_['N-1'])+1):
            self.gconst = arrset(self.gconst,ig_['2.7-'+str(I)],float(1.0))
        self.gconst = arrset(self.gconst,ig_['2.8'],float(v_['Q']))
        self.gconst  = (
              arrset(self.gconst,ig_['2.9-'+str(int(v_['K']))],float(v_['SMALLHF'])))
        self.gconst  = (
              arrset(self.gconst,ig_['2.9-'+str(int(v_['K+']))],float(v_['BIGHF'])))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.x0[ix_['X0,1']] = float(0.09203)
        self.x0[ix_['X0,2']] = float(0.908)
        self.x0[ix_['X1,1']] = float(0.1819)
        self.x0[ix_['X1,2']] = float(0.8181)
        self.x0[ix_['X2,1']] = float(0.284)
        self.x0[ix_['X2,2']] = float(0.716)
        self.x0[ix_['X3,1']] = float(0.3051)
        self.x0[ix_['X3,2']] = float(0.6949)
        self.x0[ix_['X4,1']] = float(0.3566)
        self.x0[ix_['X4,2']] = float(0.6434)
        self.x0[ix_['X5,1']] = float(0.468)
        self.x0[ix_['X5,2']] = float(0.532)
        self.x0[ix_['X6,1']] = float(0.6579)
        self.x0[ix_['X6,2']] = float(0.3421)
        self.x0[ix_['X7,1']] = float(0.8763)
        self.x0[ix_['X7,2']] = float(0.1237)
        self.x0[ix_['T0']] = float(120.0)
        self.x0[ix_['T1']] = float(110.0)
        self.x0[ix_['T2']] = float(100.0)
        self.x0[ix_['T3']] = float(88.0)
        self.x0[ix_['T4']] = float(86.0)
        self.x0[ix_['T5']] = float(84.0)
        self.x0[ix_['T6']] = float(80.0)
        self.x0[ix_['T7']] = float(76.0)
        self.x0[ix_['V0']] = float(886.37)
        self.x0[ix_['V1']] = float(910.01)
        self.x0[ix_['V2']] = float(922.52)
        self.x0[ix_['V3']] = float(926.46)
        self.x0[ix_['V4']] = float(935.56)
        self.x0[ix_['V5']] = float(952.83)
        self.x0[ix_['V6']] = float(975.73)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'en2PROD', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftp = []
        elftp = loaset(elftp,it,0,'P1')
        elftp = loaset(elftp,it,1,'P2')
        [it,iet_,_] = s2mpj_ii( 'ePOLY1PRD', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftp = loaset(elftp,it,0,'P1')
        elftp = loaset(elftp,it,1,'P6')
        elftp = loaset(elftp,it,2,'P7')
        elftp = loaset(elftp,it,3,'P8')
        [it,iet_,_] = s2mpj_ii( 'ePOLY2PRD', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftp = loaset(elftp,it,0,'P1')
        elftp = loaset(elftp,it,1,'P2')
        elftp = loaset(elftp,it,2,'P6')
        elftp = loaset(elftp,it,3,'P7')
        elftp = loaset(elftp,it,4,'P8')
        [it,iet_,_] = s2mpj_ii( 'eEXP2PROD', iet_)
        elftv = loaset(elftv,it,0,'V2')
        elftv = loaset(elftv,it,1,'V3')
        elftp = loaset(elftp,it,0,'P1')
        elftp = loaset(elftp,it,1,'P2')
        elftp = loaset(elftp,it,2,'P3')
        elftp = loaset(elftp,it,3,'P4')
        elftp = loaset(elftp,it,4,'P5')
        [it,iet_,_] = s2mpj_ii( 'eEXP3PROD', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftp = loaset(elftp,it,0,'P1')
        elftp = loaset(elftp,it,1,'P2')
        elftp = loaset(elftp,it,2,'P3')
        elftp = loaset(elftp,it,3,'P4')
        elftp = loaset(elftp,it,4,'P5')
        [it,iet_,_] = s2mpj_ii( 'eEXP4PROD', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftp = loaset(elftp,it,0,'P1')
        elftp = loaset(elftp,it,1,'P2')
        elftp = loaset(elftp,it,2,'P3')
        elftp = loaset(elftp,it,3,'P4')
        elftp = loaset(elftp,it,4,'P5')
        elftp = loaset(elftp,it,5,'P6')
        elftp = loaset(elftp,it,6,'P7')
        elftp = loaset(elftp,it,7,'P8')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        v_['-D'] = -1.0*v_['D']
        for J in range(int(v_['1']),int(v_['M'])+1):
            ename = 'E11-'+str(J)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'en2PROD')
            ielftype = arrset(ielftype,ie,iet_["en2PROD"])
            vname = 'X'+str(int(v_['1']))+','+str(J)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'V'+str(int(v_['0']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='P1')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(-1.0))
            posep = np.where(elftp[ielftype[ie]]=='P2')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['B']))
            ename = 'E12-'+str(J)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eEXP3PROD')
            ielftype = arrset(ielftype,ie,iet_["eEXP3PROD"])
            vname = 'V'+str(int(v_['0']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['0']))+','+str(J)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'T'+str(int(v_['0']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='P1')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
            posep = np.where(elftp[ielftype[ie]]=='P2')[0]
            self.elpar  = (
                  loaset(self.elpar,ie,posep[0],float(v_['INVPI'+str(int(v_['0']))])))
            posep = np.where(elftp[ielftype[ie]]=='P3')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['A'+str(J)]))
            posep = np.where(elftp[ielftype[ie]]=='P4')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['B'+str(J)]))
            posep = np.where(elftp[ielftype[ie]]=='P5')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['C'+str(J)]))
            for I in range(int(v_['1']),int(v_['N-2'])+1):
                v_['I-1'] = -1+I
                v_['I+1'] = 1+I
                ename = 'E21-'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'en2PROD')
                ielftype = arrset(ielftype,ie,iet_["en2PROD"])
                vname = 'X'+str(int(v_['I+1']))+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'V'+str(I)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                posep = np.where(elftp[ielftype[ie]]=='P1')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(-1.0))
                ename = 'E22-'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eEXP3PROD')
                ielftype = arrset(ielftype,ie,iet_["eEXP3PROD"])
                vname = 'V'+str(int(v_['I-1']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'X'+str(int(v_['I-1']))+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'T'+str(int(v_['I-1']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='V3')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                posep = np.where(elftp[ielftype[ie]]=='P1')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(-1.0))
                posep = np.where(elftp[ielftype[ie]]=='P2')[0]
                self.elpar  = (
                      loaset(self.elpar,ie,posep[0],float(v_['INVPI'+str(int(v_['I-1']))])))
                posep = np.where(elftp[ielftype[ie]]=='P3')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['A'+str(J)]))
                posep = np.where(elftp[ielftype[ie]]=='P4')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['B'+str(J)]))
                posep = np.where(elftp[ielftype[ie]]=='P5')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['C'+str(J)]))
                ename = 'E23-'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'en2PROD')
                ielftype = arrset(ielftype,ie,iet_["en2PROD"])
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'V'+str(int(v_['I-1']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                posep = np.where(elftp[ielftype[ie]]=='P1')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
                ename = 'E24-'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eEXP3PROD')
                ielftype = arrset(ielftype,ie,iet_["eEXP3PROD"])
                vname = 'V'+str(I)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'T'+str(I)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='V3')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                posep = np.where(elftp[ielftype[ie]]=='P1')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
                posep = np.where(elftp[ielftype[ie]]=='P2')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['INVPI'+str(I)]))
                posep = np.where(elftp[ielftype[ie]]=='P3')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['A'+str(J)]))
                posep = np.where(elftp[ielftype[ie]]=='P4')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['B'+str(J)]))
                posep = np.where(elftp[ielftype[ie]]=='P5')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['C'+str(J)]))
            for I in range(int(v_['1']),int(v_['K-'])+1):
                ename = 'E21-'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                posep = np.where(elftp[ielftype[ie]]=='P2')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['B']))
                ename = 'E23-'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                posep = np.where(elftp[ielftype[ie]]=='P2')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['B']))
            ename = 'E21-'+str(int(v_['K']))+','+str(J)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            posep = np.where(elftp[ielftype[ie]]=='P2')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['-D']))
            ename = 'E23-'+str(int(v_['K']))+','+str(J)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            posep = np.where(elftp[ielftype[ie]]=='P2')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['B']))
            for I in range(int(v_['K+']),int(v_['N-2'])+1):
                ename = 'E21-'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                posep = np.where(elftp[ielftype[ie]]=='P2')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['-D']))
                ename = 'E23-'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                posep = np.where(elftp[ielftype[ie]]=='P2')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['-D']))
            ename = 'E31-'+str(J)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eEXP2PROD')
            ielftype = arrset(ielftype,ie,iet_["eEXP2PROD"])
            vname = 'X'+str(int(v_['N-2']))+','+str(J)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'T'+str(int(v_['N-2']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='P1')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
            posep = np.where(elftp[ielftype[ie]]=='P2')[0]
            self.elpar  = (
                  loaset(self.elpar,ie,posep[0],float(v_['INVPI'+str(int(v_['N-2']))])))
            posep = np.where(elftp[ielftype[ie]]=='P3')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['A'+str(J)]))
            posep = np.where(elftp[ielftype[ie]]=='P4')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['B'+str(J)]))
            posep = np.where(elftp[ielftype[ie]]=='P5')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['C'+str(J)]))
        for J in range(int(v_['1']),int(v_['M'])+1):
            for I in range(int(v_['0']),int(v_['N-1'])+1):
                ename = 'E71-'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eEXP2PROD')
                ielftype = arrset(ielftype,ie,iet_["eEXP2PROD"])
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'T'+str(I)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='V3')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                posep = np.where(elftp[ielftype[ie]]=='P1')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
                posep = np.where(elftp[ielftype[ie]]=='P2')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['INVPI'+str(I)]))
                posep = np.where(elftp[ielftype[ie]]=='P3')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['A'+str(J)]))
                posep = np.where(elftp[ielftype[ie]]=='P4')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['B'+str(J)]))
                posep = np.where(elftp[ielftype[ie]]=='P5')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['C'+str(J)]))
        for J in range(int(v_['1']),int(v_['M'])+1):
            ename = 'E81-'+str(J)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eEXP4PROD')
            ielftype = arrset(ielftype,ie,iet_["eEXP4PROD"])
            vname = 'V'+str(int(v_['0']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['0']))+','+str(J)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'T'+str(int(v_['0']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='P1')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
            posep = np.where(elftp[ielftype[ie]]=='P2')[0]
            self.elpar  = (
                  loaset(self.elpar,ie,posep[0],float(v_['INVPI'+str(int(v_['0']))])))
            posep = np.where(elftp[ielftype[ie]]=='P3')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['A'+str(J)]))
            posep = np.where(elftp[ielftype[ie]]=='P4')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['B'+str(J)]))
            posep = np.where(elftp[ielftype[ie]]=='P5')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['C'+str(J)]))
            posep = np.where(elftp[ielftype[ie]]=='P6')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['BE'+str(J)]))
            posep = np.where(elftp[ielftype[ie]]=='P7')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['BEp'+str(J)]))
            posep = np.where(elftp[ielftype[ie]]=='P8')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['BEpp'+str(J)]))
            ename = 'E82-'+str(J)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'ePOLY1PRD')
            ielftype = arrset(ielftype,ie,iet_["ePOLY1PRD"])
            vname = 'X'+str(int(v_['0']))+','+str(J)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'T'+str(int(v_['0']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='P1')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['B']))
            posep = np.where(elftp[ielftype[ie]]=='P6')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['AL'+str(J)]))
            posep = np.where(elftp[ielftype[ie]]=='P7')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['ALp'+str(J)]))
            posep = np.where(elftp[ielftype[ie]]=='P8')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['ALpp'+str(J)]))
            ename = 'E83-'+str(J)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'ePOLY2PRD')
            ielftype = arrset(ielftype,ie,iet_["ePOLY2PRD"])
            vname = 'X'+str(int(v_['1']))+','+str(J)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'V'+str(int(v_['0']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'T'+str(int(v_['1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='P1')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(-1.0))
            posep = np.where(elftp[ielftype[ie]]=='P2')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['B']))
            posep = np.where(elftp[ielftype[ie]]=='P6')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['AL'+str(J)]))
            posep = np.where(elftp[ielftype[ie]]=='P7')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['ALp'+str(J)]))
            posep = np.where(elftp[ielftype[ie]]=='P8')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['ALpp'+str(J)]))
            for I in range(int(v_['1']),int(v_['N-2'])+1):
                v_['I-1'] = -1+I
                v_['I+1'] = 1+I
                ename = 'E91-'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eEXP4PROD')
                ielftype = arrset(ielftype,ie,iet_["eEXP4PROD"])
                vname = 'V'+str(I)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'T'+str(I)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='V3')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                posep = np.where(elftp[ielftype[ie]]=='P1')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
                posep = np.where(elftp[ielftype[ie]]=='P2')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['INVPI'+str(I)]))
                posep = np.where(elftp[ielftype[ie]]=='P3')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['A'+str(J)]))
                posep = np.where(elftp[ielftype[ie]]=='P4')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['B'+str(J)]))
                posep = np.where(elftp[ielftype[ie]]=='P5')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['C'+str(J)]))
                posep = np.where(elftp[ielftype[ie]]=='P6')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['BE'+str(J)]))
                posep = np.where(elftp[ielftype[ie]]=='P7')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['BEp'+str(J)]))
                posep = np.where(elftp[ielftype[ie]]=='P8')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['BEpp'+str(J)]))
                ename = 'E92-'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'ePOLY2PRD')
                ielftype = arrset(ielftype,ie,iet_["ePOLY2PRD"])
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'V'+str(int(v_['I-1']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'T'+str(I)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='V3')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                posep = np.where(elftp[ielftype[ie]]=='P1')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
                posep = np.where(elftp[ielftype[ie]]=='P6')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['AL'+str(J)]))
                posep = np.where(elftp[ielftype[ie]]=='P7')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['ALp'+str(J)]))
                posep = np.where(elftp[ielftype[ie]]=='P8')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['ALpp'+str(J)]))
                ename = 'E93-'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eEXP4PROD')
                ielftype = arrset(ielftype,ie,iet_["eEXP4PROD"])
                vname = 'V'+str(int(v_['I-1']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'X'+str(int(v_['I-1']))+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'T'+str(int(v_['I-1']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='V3')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                posep = np.where(elftp[ielftype[ie]]=='P1')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(-1.0))
                posep = np.where(elftp[ielftype[ie]]=='P2')[0]
                self.elpar  = (
                      loaset(self.elpar,ie,posep[0],float(v_['INVPI'+str(int(v_['I-1']))])))
                posep = np.where(elftp[ielftype[ie]]=='P3')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['A'+str(J)]))
                posep = np.where(elftp[ielftype[ie]]=='P4')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['B'+str(J)]))
                posep = np.where(elftp[ielftype[ie]]=='P5')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['C'+str(J)]))
                posep = np.where(elftp[ielftype[ie]]=='P6')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['BE'+str(J)]))
                posep = np.where(elftp[ielftype[ie]]=='P7')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['BEp'+str(J)]))
                posep = np.where(elftp[ielftype[ie]]=='P8')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['BEpp'+str(J)]))
                ename = 'E94-'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'ePOLY2PRD')
                ielftype = arrset(ielftype,ie,iet_["ePOLY2PRD"])
                vname = 'X'+str(int(v_['I+1']))+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'V'+str(I)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'T'+str(int(v_['I+1']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='V3')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                posep = np.where(elftp[ielftype[ie]]=='P1')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(-1.0))
                posep = np.where(elftp[ielftype[ie]]=='P6')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['AL'+str(J)]))
                posep = np.where(elftp[ielftype[ie]]=='P7')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['ALp'+str(J)]))
                posep = np.where(elftp[ielftype[ie]]=='P8')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['ALpp'+str(J)]))
            for I in range(int(v_['1']),int(v_['K-'])+1):
                ename = 'E92-'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                posep = np.where(elftp[ielftype[ie]]=='P2')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['B']))
                ename = 'E94-'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                posep = np.where(elftp[ielftype[ie]]=='P2')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['B']))
            ename = 'E92-'+str(int(v_['K']))+','+str(J)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            posep = np.where(elftp[ielftype[ie]]=='P2')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['B']))
            ename = 'E94-'+str(int(v_['K']))+','+str(J)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            posep = np.where(elftp[ielftype[ie]]=='P2')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['-D']))
            for I in range(int(v_['K+']),int(v_['N-2'])+1):
                ename = 'E92-'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                posep = np.where(elftp[ielftype[ie]]=='P2')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['-D']))
                ename = 'E94-'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                posep = np.where(elftp[ielftype[ie]]=='P2')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['-D']))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gL2',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for ig in range(0,ngrp):
            self.grftype = arrset(self.grftype,ig,'gL2')
        for J in range(int(v_['1']),int(v_['M'])+1):
            ig = ig_['2.1-'+str(J)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['E11-'+str(J)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = posel+1
            self.grelt = loaset(self.grelt,ig,posel,ie_['E12-'+str(J)])
            self.grelw = loaset(self.grelw,ig,posel, 1.)
            ig = ig_['2.3-'+str(J)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['E31-'+str(J)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
            ig = ig_['2.8']
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['E81-'+str(J)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = posel+1
            self.grelt = loaset(self.grelt,ig,posel,ie_['E82-'+str(J)])
            self.grelw = loaset(self.grelw,ig,posel, 1.)
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['E83-'+str(J)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
            for I in range(int(v_['1']),int(v_['N-2'])+1):
                ig = ig_['2.2-'+str(I)+','+str(J)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['E21-'+str(I)+','+str(J)])
                self.grelw = loaset(self.grelw,ig,posel,1.)
                posel = posel+1
                self.grelt = loaset(self.grelt,ig,posel,ie_['E22-'+str(I)+','+str(J)])
                self.grelw = loaset(self.grelw,ig,posel, 1.)
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['E23-'+str(I)+','+str(J)])
                self.grelw = loaset(self.grelw,ig,posel,1.)
                posel = posel+1
                self.grelt = loaset(self.grelt,ig,posel,ie_['E24-'+str(I)+','+str(J)])
                self.grelw = loaset(self.grelw,ig,posel, 1.)
                ig = ig_['2.9-'+str(I)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['E91-'+str(I)+','+str(J)])
                self.grelw = loaset(self.grelw,ig,posel,1.)
                posel = posel+1
                self.grelt = loaset(self.grelt,ig,posel,ie_['E92-'+str(I)+','+str(J)])
                self.grelw = loaset(self.grelw,ig,posel, 1.)
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['E93-'+str(I)+','+str(J)])
                self.grelw = loaset(self.grelw,ig,posel,1.)
                posel = posel+1
                self.grelt = loaset(self.grelt,ig,posel,ie_['E94-'+str(I)+','+str(J)])
                self.grelw = loaset(self.grelw,ig,posel, 1.)
            for I in range(int(v_['0']),int(v_['N-1'])+1):
                ig = ig_['2.7-'+str(I)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['E71-'+str(I)+','+str(J)])
                self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        self.objlower = 0.0
#    Solution
# LO SOLTN               0.0
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass   = "C-CSUR2-MN-31-0"
        self.objderlvl = 2

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def en2PROD(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = self.elpar[iel_][0]*EV_[0]*(EV_[1]+self.elpar[iel_][1])
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = self.elpar[iel_][0]*(EV_[1]+self.elpar[iel_][1])
            g_[1] = self.elpar[iel_][0]*EV_[0]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = self.elpar[iel_][0]
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ePOLY1PRD(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        POLY  = (
              self.elpar[iel_][1]+self.elpar[iel_][2]*EV_[1]+self.elpar[iel_][3]*EV_[1]*EV_[1])
        DPOLY = self.elpar[iel_][2]+2.0*self.elpar[iel_][3]*EV_[1]
        f_   = self.elpar[iel_][0]*EV_[0]*POLY
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = self.elpar[iel_][0]*POLY
            g_[1] = self.elpar[iel_][0]*EV_[0]*DPOLY
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = self.elpar[iel_][0]*DPOLY
                H_[1,0] = H_[0,1]
                H_[1,1] = self.elpar[iel_][0]*EV_[0]*2.0e+0*self.elpar[iel_][3]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ePOLY2PRD(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        POLY  = (
              self.elpar[iel_][2]+self.elpar[iel_][3]*EV_[2]+self.elpar[iel_][4]*EV_[2]*EV_[2])
        DPOLY = self.elpar[iel_][3]+2.0*self.elpar[iel_][4]*EV_[2]
        f_   = self.elpar[iel_][0]*EV_[0]*(self.elpar[iel_][1]+EV_[1])*POLY
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = self.elpar[iel_][0]*(self.elpar[iel_][1]+EV_[1])*POLY
            g_[1] = self.elpar[iel_][0]*EV_[0]*POLY
            g_[2] = self.elpar[iel_][0]*EV_[0]*(self.elpar[iel_][1]+EV_[1])*DPOLY
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = self.elpar[iel_][0]*POLY
                H_[1,0] = H_[0,1]
                H_[0,2] = self.elpar[iel_][0]*(self.elpar[iel_][1]+EV_[1])*DPOLY
                H_[2,0] = H_[0,2]
                H_[1,2] = self.elpar[iel_][0]*EV_[0]*DPOLY
                H_[2,1] = H_[1,2]
                H_[2,2]  = (
                      self.elpar[iel_][0]*EV_[0]*(self.elpar[iel_][1]+EV_[1])*2.0e+0*self.elpar[iel_][4])
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eEXP2PROD(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        EXPROD  = (
              self.elpar[iel_][0]*self.elpar[iel_][1]*np.exp(self.elpar[iel_][2]+(self.elpar[iel_][3]/(EV_[1]+self.elpar[iel_][4]))))
        F = EV_[0]*EXPROD
        f_   = F
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EXPROD
            g_[1] = -EV_[0]*EXPROD*self.elpar[iel_][3]/(self.elpar[iel_][4]+EV_[1])**2
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = -EXPROD*self.elpar[iel_][3]/(self.elpar[iel_][4]+EV_[1])**2
                H_[1,0] = H_[0,1]
                H_[1,1] = (F*(self.elpar[iel_][3]/(self.elpar[iel_][4]+EV_[1])**2)**2+
                     2.0e+0*F*self.elpar[iel_][3]/(self.elpar[iel_][4]+EV_[1])**3)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eEXP3PROD(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        EXPROD  = (
              self.elpar[iel_][0]*self.elpar[iel_][1]*np.exp(self.elpar[iel_][2]+(self.elpar[iel_][3]/(EV_[2]+self.elpar[iel_][4]))))
        F = EV_[0]*EV_[1]*EXPROD
        TERM = -self.elpar[iel_][3]/(self.elpar[iel_][4]+EV_[2])**2
        f_   = F
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]*EXPROD
            g_[1] = EV_[0]*EXPROD
            g_[2] = F*TERM
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = EXPROD
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1]*EXPROD*TERM
                H_[2,0] = H_[0,2]
                H_[1,2] = EV_[0]*EXPROD*TERM
                H_[2,1] = H_[1,2]
                H_[2,2]  = (
                      F*(TERM*TERM+2.0e+0*self.elpar[iel_][3]/(self.elpar[iel_][4]+EV_[2])**3))
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eEXP4PROD(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        EXPROD  = (
              self.elpar[iel_][0]*self.elpar[iel_][1]*np.exp(self.elpar[iel_][2]+(self.elpar[iel_][3]/(EV_[2]+self.elpar[iel_][4]))))
        F = EV_[0]*EV_[1]*EXPROD
        POLY  = (
              self.elpar[iel_][5]+self.elpar[iel_][6]*EV_[2]+self.elpar[iel_][7]*EV_[2]*EV_[2])
        DPOLY = self.elpar[iel_][6]+2.0*self.elpar[iel_][7]*EV_[2]
        TERM = DPOLY-POLY*self.elpar[iel_][3]/(self.elpar[iel_][4]+EV_[2])**2
        f_   = F*POLY
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]*EXPROD*POLY
            g_[1] = EV_[0]*EXPROD*POLY
            g_[2] = F*TERM
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = EXPROD*POLY
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1]*EXPROD*TERM
                H_[2,0] = H_[0,2]
                H_[1,2] = EV_[0]*EXPROD*TERM
                H_[2,1] = H_[1,2]
                H_[2,2]  = (
                      F*(-(self.elpar[iel_][3]/(self.elpar[iel_][4]+EV_[2])**2)*TERM+2.0*self.elpar[iel_][7]-DPOLY*self.elpar[iel_][3]/(self.elpar[iel_][4]+EV_[2])**2+2.0e+0*POLY*self.elpar[iel_][3]/(self.elpar[iel_][4]+EV_[2])**3))
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gL2(self,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= GVAR_*GVAR_
        if nargout>1:
            g_ = GVAR_+GVAR_
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

