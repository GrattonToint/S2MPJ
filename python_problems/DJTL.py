from s2mpjlib import *
class  DJTL(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DJTL
#    *********
# 
#    Source: modified version of problem 19 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
#    that is meant to simulate the Lagrangian barrier objective function
#    for particular values of the shifts and multipliers
# 
#    SIF input: A.R. Conn August 1993
# 
#    classification = "C-COUR2-AN-2-0"
# 
#    Define multipliers and shifts
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'DJTL'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['LL1'] = 1.0
        v_['LL2'] = 1.0
        v_['LL3'] = 1.0
        v_['LL4'] = 1.0
        v_['SL1'] = 1.0
        v_['SL2'] = 1.0
        v_['SL3'] = 1.0
        v_['SL4'] = 1.0
        v_['LU1'] = 1.0
        v_['LU2'] = 1.0
        v_['LU3'] = 1.0
        v_['LU4'] = 1.0
        v_['SU1'] = 1.0
        v_['SU2'] = 1.0
        v_['SU3'] = 1.0
        v_['SU4'] = 1.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [iv,ix_,_] = s2mpj_ii('X1',ix_)
        self.xnames=arrset(self.xnames,iv,'X1')
        [iv,ix_,_] = s2mpj_ii('X2',ix_)
        self.xnames=arrset(self.xnames,iv,'X2')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('CONU1',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('CONL1',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('CONU2',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('CONL2',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('BNDU1',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X1']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('BNDL1',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X1']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('BNDU2',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('BNDL2',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2']])
        valA = np.append(valA,float(1.0))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        selfnob      = ngrp
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        self.gconst = arrset(self.gconst,ig_['CONU1'],float(-200.0))
        self.gconst = arrset(self.gconst,ig_['CONL1'],float(100.0))
        self.gconst = arrset(self.gconst,ig_['CONU2'],float(0.0))
        self.gconst = arrset(self.gconst,ig_['CONL2'],float(-82.81))
        self.gconst = arrset(self.gconst,ig_['BNDU1'],float(-100.0))
        self.gconst = arrset(self.gconst,ig_['BNDL1'],float(13.0))
        self.gconst = arrset(self.gconst,ig_['BNDU2'],float(-100.0))
        self.gconst = arrset(self.gconst,ig_['BNDL2'],float(0.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.x0[ix_['X1']] = float(15.0)
        self.x0[ix_['X2']] = float(6.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eCBm10', iet_)
        elftv = loaset(elftv,it,0,'V1')
        [it,iet_,_] = s2mpj_ii( 'eCBm20', iet_)
        elftv = loaset(elftv,it,0,'V1')
        [it,iet_,_] = s2mpj_ii( 'eSQm5', iet_)
        elftv = loaset(elftv,it,0,'V1')
        [it,iet_,_] = s2mpj_ii( 'eSQm6', iet_)
        elftv = loaset(elftv,it,0,'V1')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        ename = 'E1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCBm10')
        ielftype = arrset(ielftype,ie,iet_["eCBm10"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCBm20')
        ielftype = arrset(ielftype,ie,iet_["eCBm20"])
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQm5')
        ielftype = arrset(ielftype,ie,iet_["eSQm5"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQm5')
        ielftype = arrset(ielftype,ie,iet_["eSQm5"])
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E5'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQm6')
        ielftype = arrset(ielftype,ie,iet_["eSQm6"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gLOG',igt_)
        [it,igt_,_] = s2mpj_ii('gLOG',igt_)
        grftp = []
        grftp = loaset(grftp,it,0,'P1')
        [it,igt_,_] = s2mpj_ii('gLOG',igt_)
        grftp = loaset(grftp,it,1,'P2')
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        self.grpar   = []
        ig = ig_['OBJ']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E1'])
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E2'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        ig = ig_['CONL1']
        self.grftype = arrset(self.grftype,ig,'gLOG')
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E3'])
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E4'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posgp = np.where(grftp[igt_[self.grftype[ig]]]=='P1')[0]
        self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['SL1']))
        posgp = np.where(grftp[igt_[self.grftype[ig]]]=='P2')[0]
        self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['LL1']))
        ig = ig_['CONU1']
        self.grftype = arrset(self.grftype,ig,'gLOG')
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E3'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E4'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posgp = np.where(grftp[igt_[self.grftype[ig]]]=='P1')[0]
        self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['SU1']))
        posgp = np.where(grftp[igt_[self.grftype[ig]]]=='P2')[0]
        self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['LU1']))
        ig = ig_['CONL2']
        self.grftype = arrset(self.grftype,ig,'gLOG')
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E4'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E5'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posgp = np.where(grftp[igt_[self.grftype[ig]]]=='P1')[0]
        self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['SL2']))
        posgp = np.where(grftp[igt_[self.grftype[ig]]]=='P2')[0]
        self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['LL2']))
        ig = ig_['CONU2']
        self.grftype = arrset(self.grftype,ig,'gLOG')
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E4'])
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E5'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posgp = np.where(grftp[igt_[self.grftype[ig]]]=='P1')[0]
        self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['SU2']))
        posgp = np.where(grftp[igt_[self.grftype[ig]]]=='P2')[0]
        self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['LU2']))
        ig = ig_['BNDL1']
        self.grftype = arrset(self.grftype,ig,'gLOG')
        posgp = np.where(grftp[igt_[self.grftype[ig]]]=='P1')[0]
        self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['SL3']))
        posgp = np.where(grftp[igt_[self.grftype[ig]]]=='P2')[0]
        self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['LL3']))
        ig = ig_['BNDU1']
        self.grftype = arrset(self.grftype,ig,'gLOG')
        posgp = np.where(grftp[igt_[self.grftype[ig]]]=='P1')[0]
        self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['SU3']))
        posgp = np.where(grftp[igt_[self.grftype[ig]]]=='P2')[0]
        self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['LU3']))
        ig = ig_['BNDL2']
        self.grftype = arrset(self.grftype,ig,'gLOG')
        posgp = np.where(grftp[igt_[self.grftype[ig]]]=='P1')[0]
        self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['SL4']))
        posgp = np.where(grftp[igt_[self.grftype[ig]]]=='P2')[0]
        self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['LL4']))
        ig = ig_['BNDU2']
        self.grftype = arrset(self.grftype,ig,'gLOG')
        posgp = np.where(grftp[igt_[self.grftype[ig]]]=='P1')[0]
        self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['SU4']))
        posgp = np.where(grftp[igt_[self.grftype[ig]]]=='P2')[0]
        self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['LU4']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               -8951.54472
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass   = "C-COUR2-AN-2-0"
        self.objderlvl = 2

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eCBm10(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        DIF = EV_[0,0]-10.0
        f_   = DIF**3
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 3.0*DIF*DIF
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 6.0*DIF
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eCBm20(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        DIF = EV_[0,0]-20.0
        f_   = DIF**3
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 3.0*DIF*DIF
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 6.0*DIF
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eSQm5(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        DIF = EV_[0,0]-5.0
        f_   = DIF**2
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*DIF
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eSQm6(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        DIF = EV_[0,0]-6.0
        f_   = DIF**2
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*DIF
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gLOG(self,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        APP1 = GVAR_+self.grpar[igr_][0]
        P1P2 = self.grpar[igr_][0]*self.grpar[igr_][1]
        ARG0 = APP1<=0.0
        BIG = 1.0000e+10
        if ARG0!=0:
            FF = BIG*GVAR_**2
        if ARG0==0:
            FF = -P1P2*np.log(APP1)
        if ARG0!=0:
            GG = 2.0*BIG*GVAR_
        if ARG0==0:
            GG = -P1P2/APP1
        if ARG0!=0:
            HH = 2.0*BIG
        if ARG0==0:
            HH = P1P2/APP1**2
        f_= FF
        if nargout>1:
            g_ = GG
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = HH
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

