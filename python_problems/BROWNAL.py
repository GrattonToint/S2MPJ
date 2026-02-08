from s2mpjlib import *
class  BROWNAL(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : BROWNAL
#    *********
#    Brown almost linear least squares problem.
#    This problem is a sum of n least-squares groups, the last one of
#    which has a nonlinear element.
#    It Hessian matrix is dense.
# 
#    Source: Problem 27 in
#    J.J. More', B.S. Garbow and K.E. Hillstrom,
#    "Testing Unconstrained Optimization Software",
#    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
# 
#    See also Buckley#79
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "C-CSUR2-AN-V-0"
# 
#    N is the number of free variables (variable).
# 
#           Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER     original value
# IE N                   100            $-PARAMETER
# IE N                   200            $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'BROWNAL'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(10);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
# IE N                   1000           $-PARAMETER
        v_['1'] = 1
        v_['N-1'] = -1+v_['N']
        v_['N+1'] = 1+v_['N']
        v_['RN+1'] = float(v_['N+1'])
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
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            for J in range(int(v_['1']),int(v_['I-1'])+1):
                [ig,ig_,_] = s2mpj_ii('G'+str(I),ig_)
                gtype = arrset(gtype,ig,'<>')
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(J)]])
                valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('G'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(I)]])
            valA = np.append(valA,float(2.0))
            for J in range(int(v_['I+1']),int(v_['N'])+1):
                [ig,ig_,_] = s2mpj_ii('G'+str(I),ig_)
                gtype = arrset(gtype,ig,'<>')
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(J)]])
                valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('G'+str(int(v_['N'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        selfnob      = ngrp
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            self.gconst = arrset(self.gconst,ig_['G'+str(I)],float(v_['RN+1']))
        self.gconst = arrset(self.gconst,ig_['G'+str(int(v_['N']))],float(1.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        for I in range(int(v_['1']),int(v_['N'])+1):
            self.x0[ix_['X'+str(I)]] = float(0.5)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'ePROD', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftv = loaset(elftv,it,3,'V4')
        elftv = loaset(elftv,it,4,'V5')
        elftv = loaset(elftv,it,5,'V6')
        elftv = loaset(elftv,it,6,'V7')
        elftv = loaset(elftv,it,7,'V8')
        elftv = loaset(elftv,it,8,'V9')
        elftv = loaset(elftv,it,9,'V10')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        ename = 'E'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePROD')
        ielftype = arrset(ielftype,ie,iet_["ePROD"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V4')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V5')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V6')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V7')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X8'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V8')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X9'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V9')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X10'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V10')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
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
        ig = ig_['G'+str(int(v_['N']))]
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E'])
        self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 0.0
#    Solution
# LO SOLTN               0.0
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass   = "C-CSUR2-AN-V-0"
        self.objderlvl = 2

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def ePROD(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        V12 = EV_[0,0]*EV_[1,0]
        V34 = EV_[2,0]*EV_[3,0]
        V56 = EV_[4,0]*EV_[5,0]
        V78 = EV_[6,0]*EV_[7,0]
        V910 = EV_[8,0]*EV_[9,0]
        f_   = V12*V34*V56*V78*V910
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1,0]*V34*V56*V78*V910
            g_[1] = EV_[0,0]*V34*V56*V78*V910
            g_[2] = V12*EV_[3,0]*V56*V78*V910
            g_[3] = V12*EV_[2,0]*V56*V78*V910
            g_[4] = V12*V34*EV_[5,0]*V78*V910
            g_[5] = V12*V34*EV_[4,0]*V78*V910
            g_[6] = V12*V34*V56*EV_[7,0]*V910
            g_[7] = V12*V34*V56*EV_[6,0]*V910
            g_[8] = V12*V34*V56*V78*EV_[9,0]
            g_[9] = V12*V34*V56*V78*EV_[8,0]
            if nargout>2:
                H_ = np.zeros((10,10))
                H_[0,1] = V34*V56*V78*V910
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1,0]*EV_[3,0]*V56*V78*V910
                H_[2,0] = H_[0,2]
                H_[0,3] = EV_[1,0]*EV_[2,0]*V56*V78*V910
                H_[3,0] = H_[0,3]
                H_[0,4] = EV_[1,0]*V34*EV_[5,0]*V78*V910
                H_[4,0] = H_[0,4]
                H_[0,5] = EV_[1,0]*V34*EV_[4,0]*V78*V910
                H_[5,0] = H_[0,5]
                H_[0,6] = EV_[1,0]*V34*V56*EV_[7,0]*V910
                H_[6,0] = H_[0,6]
                H_[0,7] = EV_[1,0]*V34*V56*EV_[6,0]*V910
                H_[7,0] = H_[0,7]
                H_[0,8] = EV_[1,0]*V34*V56*V78*EV_[9,0]
                H_[8,0] = H_[0,8]
                H_[0,9] = EV_[1,0]*V34*V56*V78*EV_[8,0]
                H_[9,0] = H_[0,9]
                H_[1,2] = EV_[0,0]*EV_[3,0]*V56*V78*V910
                H_[2,1] = H_[1,2]
                H_[1,3] = EV_[0,0]*EV_[2,0]*V56*V78*V910
                H_[3,1] = H_[1,3]
                H_[1,4] = EV_[0,0]*V34*EV_[5,0]*V78*V910
                H_[4,1] = H_[1,4]
                H_[1,5] = EV_[0,0]*V34*EV_[4,0]*V78*V910
                H_[5,1] = H_[1,5]
                H_[1,6] = EV_[0,0]*V34*V56*EV_[7,0]*V910
                H_[6,1] = H_[1,6]
                H_[1,7] = EV_[0,0]*V34*V56*EV_[6,0]*V910
                H_[7,1] = H_[1,7]
                H_[1,8] = EV_[0,0]*V34*V56*V78*EV_[9,0]
                H_[8,1] = H_[1,8]
                H_[1,9] = EV_[0,0]*V34*V56*V78*EV_[8,0]
                H_[9,1] = H_[1,9]
                H_[2,3] = V12*V56*V78*V910
                H_[3,2] = H_[2,3]
                H_[2,4] = V12*EV_[3,0]*EV_[5,0]*V78*V910
                H_[4,2] = H_[2,4]
                H_[2,5] = V12*EV_[3,0]*EV_[4,0]*V78*V910
                H_[5,2] = H_[2,5]
                H_[2,6] = V12*EV_[3,0]*V56*EV_[7,0]*V910
                H_[6,2] = H_[2,6]
                H_[2,7] = V12*EV_[3,0]*V56*EV_[6,0]*V910
                H_[7,2] = H_[2,7]
                H_[2,8] = V12*EV_[3,0]*V56*V78*EV_[9,0]
                H_[8,2] = H_[2,8]
                H_[2,9] = V12*EV_[3,0]*V56*V78*EV_[8,0]
                H_[9,2] = H_[2,9]
                H_[3,4] = V12*EV_[2,0]*EV_[5,0]*V78*V910
                H_[4,3] = H_[3,4]
                H_[3,5] = V12*EV_[2,0]*EV_[4,0]*V78*V910
                H_[5,3] = H_[3,5]
                H_[3,6] = V12*EV_[2,0]*V56*EV_[7,0]*V910
                H_[6,3] = H_[3,6]
                H_[3,7] = V12*EV_[2,0]*V56*EV_[6,0]*V910
                H_[7,3] = H_[3,7]
                H_[3,8] = V12*EV_[2,0]*V56*V78*EV_[9,0]
                H_[8,3] = H_[3,8]
                H_[3,9] = V12*EV_[2,0]*V56*V78*EV_[8,0]
                H_[9,3] = H_[3,9]
                H_[4,5] = V12*V34*V78*V910
                H_[5,4] = H_[4,5]
                H_[4,6] = V12*V34*EV_[5,0]*EV_[7,0]*V910
                H_[6,4] = H_[4,6]
                H_[4,7] = V12*V34*EV_[5,0]*EV_[6,0]*V910
                H_[7,4] = H_[4,7]
                H_[4,8] = V12*V34*EV_[5,0]*V78*EV_[9,0]
                H_[8,4] = H_[4,8]
                H_[4,9] = V12*V34*EV_[5,0]*V78*EV_[8,0]
                H_[9,4] = H_[4,9]
                H_[5,6] = V12*V34*EV_[4,0]*EV_[7,0]*V910
                H_[6,5] = H_[5,6]
                H_[5,7] = V12*V34*EV_[4,0]*EV_[6,0]*V910
                H_[7,5] = H_[5,7]
                H_[5,8] = V12*V34*EV_[4,0]*V78*EV_[9,0]
                H_[8,5] = H_[5,8]
                H_[5,9] = V12*V34*EV_[4,0]*V78*EV_[8,0]
                H_[9,5] = H_[5,9]
                H_[6,7] = V12*V34*V56*V910
                H_[7,6] = H_[6,7]
                H_[6,8] = V12*V34*V56*EV_[7,0]*EV_[9,0]
                H_[8,6] = H_[6,8]
                H_[6,9] = V12*V34*V56*EV_[7,0]*EV_[8,0]
                H_[9,6] = H_[6,9]
                H_[7,8] = V12*V34*V56*EV_[6,0]*EV_[9,0]
                H_[8,7] = H_[7,8]
                H_[7,9] = V12*V34*V56*EV_[6,0]*EV_[8,0]
                H_[9,7] = H_[7,9]
                H_[8,9] = V12*V34*V56*V78
                H_[9,8] = H_[8,9]
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

