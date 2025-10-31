from s2mpjlib import *
class  OSBORNEB(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : OSBORNEB
#    *********
# 
#    Osborne second problem in 11 variables.
# 
#    This function  is a nonlinear least squares with 65 groups.  Each
#    group has 4 nonlinear elements.
# 
#    Source:  Problem 19 in
#    J.J. More', B.S. Garbow and K.E. Hillstrom,
#    "Testing Unconstrained Optimization Software",
#    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
# 
#    See also Buckley#32 (p.78).
# 
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "C-CSUR2-MN-11-0"
# 
#    Number of groups
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'OSBORNEB'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['M'] = 65
        v_['N'] = 11
        v_['1'] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['1']),int(v_['M'])+1):
            [ig,ig_,_] = s2mpj_ii('G'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        self.gconst = arrset(self.gconst,ig_['G1'],float(1.366))
        self.gconst = arrset(self.gconst,ig_['G2'],float(1.191))
        self.gconst = arrset(self.gconst,ig_['G3'],float(1.112))
        self.gconst = arrset(self.gconst,ig_['G4'],float(1.013))
        self.gconst = arrset(self.gconst,ig_['G5'],float(0.991))
        self.gconst = arrset(self.gconst,ig_['G6'],float(0.885))
        self.gconst = arrset(self.gconst,ig_['G7'],float(0.831))
        self.gconst = arrset(self.gconst,ig_['G8'],float(0.847))
        self.gconst = arrset(self.gconst,ig_['G9'],float(0.786))
        self.gconst = arrset(self.gconst,ig_['G10'],float(0.725))
        self.gconst = arrset(self.gconst,ig_['G11'],float(0.746))
        self.gconst = arrset(self.gconst,ig_['G12'],float(0.679))
        self.gconst = arrset(self.gconst,ig_['G13'],float(0.608))
        self.gconst = arrset(self.gconst,ig_['G14'],float(0.655))
        self.gconst = arrset(self.gconst,ig_['G15'],float(0.616))
        self.gconst = arrset(self.gconst,ig_['G16'],float(0.606))
        self.gconst = arrset(self.gconst,ig_['G17'],float(0.602))
        self.gconst = arrset(self.gconst,ig_['G18'],float(0.626))
        self.gconst = arrset(self.gconst,ig_['G19'],float(0.651))
        self.gconst = arrset(self.gconst,ig_['G20'],float(0.724))
        self.gconst = arrset(self.gconst,ig_['G21'],float(0.649))
        self.gconst = arrset(self.gconst,ig_['G22'],float(0.649))
        self.gconst = arrset(self.gconst,ig_['G23'],float(0.694))
        self.gconst = arrset(self.gconst,ig_['G24'],float(0.644))
        self.gconst = arrset(self.gconst,ig_['G25'],float(0.624))
        self.gconst = arrset(self.gconst,ig_['G26'],float(0.661))
        self.gconst = arrset(self.gconst,ig_['G27'],float(0.612))
        self.gconst = arrset(self.gconst,ig_['G28'],float(0.558))
        self.gconst = arrset(self.gconst,ig_['G29'],float(0.533))
        self.gconst = arrset(self.gconst,ig_['G30'],float(0.495))
        self.gconst = arrset(self.gconst,ig_['G31'],float(0.500))
        self.gconst = arrset(self.gconst,ig_['G32'],float(0.423))
        self.gconst = arrset(self.gconst,ig_['G33'],float(0.395))
        self.gconst = arrset(self.gconst,ig_['G34'],float(0.375))
        self.gconst = arrset(self.gconst,ig_['G35'],float(0.372))
        self.gconst = arrset(self.gconst,ig_['G36'],float(0.391))
        self.gconst = arrset(self.gconst,ig_['G37'],float(0.396))
        self.gconst = arrset(self.gconst,ig_['G38'],float(0.405))
        self.gconst = arrset(self.gconst,ig_['G39'],float(0.428))
        self.gconst = arrset(self.gconst,ig_['G40'],float(0.429))
        self.gconst = arrset(self.gconst,ig_['G41'],float(0.523))
        self.gconst = arrset(self.gconst,ig_['G42'],float(0.562))
        self.gconst = arrset(self.gconst,ig_['G43'],float(0.607))
        self.gconst = arrset(self.gconst,ig_['G44'],float(0.653))
        self.gconst = arrset(self.gconst,ig_['G45'],float(0.672))
        self.gconst = arrset(self.gconst,ig_['G46'],float(0.708))
        self.gconst = arrset(self.gconst,ig_['G47'],float(0.633))
        self.gconst = arrset(self.gconst,ig_['G48'],float(0.668))
        self.gconst = arrset(self.gconst,ig_['G49'],float(0.645))
        self.gconst = arrset(self.gconst,ig_['G50'],float(0.632))
        self.gconst = arrset(self.gconst,ig_['G51'],float(0.591))
        self.gconst = arrset(self.gconst,ig_['G52'],float(0.559))
        self.gconst = arrset(self.gconst,ig_['G53'],float(0.597))
        self.gconst = arrset(self.gconst,ig_['G54'],float(0.625))
        self.gconst = arrset(self.gconst,ig_['G55'],float(0.739))
        self.gconst = arrset(self.gconst,ig_['G56'],float(0.710))
        self.gconst = arrset(self.gconst,ig_['G57'],float(0.729))
        self.gconst = arrset(self.gconst,ig_['G58'],float(0.720))
        self.gconst = arrset(self.gconst,ig_['G59'],float(0.636))
        self.gconst = arrset(self.gconst,ig_['G60'],float(0.581))
        self.gconst = arrset(self.gconst,ig_['G61'],float(0.428))
        self.gconst = arrset(self.gconst,ig_['G62'],float(0.292))
        self.gconst = arrset(self.gconst,ig_['G63'],float(0.162))
        self.gconst = arrset(self.gconst,ig_['G64'],float(0.098))
        self.gconst = arrset(self.gconst,ig_['G65'],float(0.054))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.x0[ix_['X1']] = float(1.3)
        self.x0[ix_['X2']] = float(0.65)
        self.x0[ix_['X3']] = float(0.65)
        self.x0[ix_['X4']] = float(0.7)
        self.x0[ix_['X5']] = float(0.6)
        self.x0[ix_['X6']] = float(3.0)
        self.x0[ix_['X7']] = float(5.0)
        self.x0[ix_['X8']] = float(7.0)
        self.x0[ix_['X9']] = float(2.0)
        self.x0[ix_['X10']] = float(4.5)
        self.x0[ix_['X11']] = float(5.5)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'ePEXP', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftp = []
        elftp = loaset(elftp,it,0,'T')
        [it,iet_,_] = s2mpj_ii( 'ePEXP3', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftp = loaset(elftp,it,0,'T3')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['1']),int(v_['M'])+1):
            v_['I-1'] = 1+I
            v_['RI-1'] = float(v_['I-1'])
            v_['TI'] = 0.1*v_['RI-1']
            ename = 'A'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'ePEXP')
            ielftype = arrset(ielftype,ie,iet_["ePEXP"])
            vname = 'X1'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X5'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='T')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['TI']))
            ename = 'B'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'ePEXP3')
            ielftype = arrset(ielftype,ie,iet_["ePEXP3"])
            vname = 'X2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X9'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X6'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='T3')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['TI']))
            ename = 'C'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'ePEXP3')
            ielftype = arrset(ielftype,ie,iet_["ePEXP3"])
            vname = 'X3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X10'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X7'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='T3')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['TI']))
            ename = 'D'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'ePEXP3')
            ielftype = arrset(ielftype,ie,iet_["ePEXP3"])
            vname = 'X4'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X11'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X8'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='T3')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['TI']))
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
        for I in range(int(v_['1']),int(v_['M'])+1):
            ig = ig_['G'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['A'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = posel+1
            self.grelt = loaset(self.grelt,ig,posel,ie_['B'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel, 1.)
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['C'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = posel+1
            self.grelt = loaset(self.grelt,ig,posel,ie_['D'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel, 1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        self.objlower = 0.0
#    Solution
# LO SOLTN               0.04013774
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass   = "C-CSUR2-MN-11-0"
        self.objderlvl = 2

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def ePEXP(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        EXPA = np.exp(-self.elpar[iel_][0]*EV_[1])
        FVAL = EV_[0]*EXPA
        f_   = FVAL
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EXPA
            g_[1] = -self.elpar[iel_][0]*FVAL
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = -self.elpar[iel_][0]*EXPA
                H_[1,0] = H_[0,1]
                H_[1,1] = self.elpar[iel_][0]*self.elpar[iel_][0]*FVAL
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ePEXP3(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        TMV2 = self.elpar[iel_][0]-EV_[1]
        TMV2SQ = TMV2*TMV2
        EXPA = np.exp(-TMV2SQ*EV_[2])
        FVAL = EV_[0]*EXPA
        A = 2.0*TMV2*EV_[2]
        f_   = FVAL
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EXPA
            g_[1] = A*FVAL
            g_[2] = -TMV2SQ*FVAL
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = A*EXPA
                H_[1,0] = H_[0,1]
                H_[0,2] = -TMV2SQ*EXPA
                H_[2,0] = H_[0,2]
                H_[1,1] = (A*A-2.0*EV_[2])*FVAL
                H_[1,2] = (2.0*TMV2-A*TMV2SQ)*FVAL
                H_[2,1] = H_[1,2]
                H_[2,2] = TMV2SQ*TMV2SQ*FVAL
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

