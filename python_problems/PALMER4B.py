from s2mpjlib import *
class  PALMER4B(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : PALMER4B
#    *********
# 
#    A nonlinear least squares problem with bounds
#    arising from chemical kinetics.
# 
#    model: H-N=C=Se TZVP + MP2
#    fitting Y to A2 X**2 + A4 X**4
#                 + B / ( C + X**2 ), B, C nonnegative.
# 
#    Source:
#    M.  Palmer, Edinburgh, private communication.
# 
#    SIF input: Nick Gould, 1990.
# 
#    classification = "C-CSBR2-RN-4-0"
# 
#    Number of data points
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'PALMER4B'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['M'] = 23
        v_['1'] = 1
        v_['X1'] = -1.658063
        v_['X2'] = -1.570796
        v_['X3'] = -1.396263
        v_['X4'] = -1.221730
        v_['X5'] = -1.047198
        v_['X6'] = -0.872665
        v_['X7'] = -0.741119
        v_['X8'] = -0.698132
        v_['X9'] = -0.523599
        v_['X10'] = -0.349066
        v_['X11'] = -0.174533
        v_['X12'] = 0.0
        v_['X13'] = 0.174533
        v_['X14'] = 0.349066
        v_['X15'] = 0.523599
        v_['X16'] = 0.698132
        v_['X17'] = 0.741119
        v_['X18'] = 0.872665
        v_['X19'] = 1.047198
        v_['X20'] = 1.221730
        v_['X21'] = 1.396263
        v_['X22'] = 1.570796
        v_['X23'] = 1.658063
        v_['Y1'] = 67.27625
        v_['Y2'] = 52.8537
        v_['Y3'] = 30.2718
        v_['Y4'] = 14.9888
        v_['Y5'] = 5.5675
        v_['Y6'] = 0.92603
        v_['Y7'] = 0.0
        v_['Y8'] = 0.085108
        v_['Y9'] = 1.867422
        v_['Y10'] = 5.014768
        v_['Y11'] = 8.263520
        v_['Y12'] = 9.8046208
        v_['Y13'] = 8.263520
        v_['Y14'] = 5.014768
        v_['Y15'] = 1.867422
        v_['Y16'] = 0.085108
        v_['Y17'] = 0.0
        v_['Y18'] = 0.92603
        v_['Y19'] = 5.5675
        v_['Y20'] = 14.9888
        v_['Y21'] = 30.2718
        v_['Y22'] = 52.8537
        v_['Y23'] = 67.27625
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [iv,ix_,_] = s2mpj_ii('A2',ix_)
        self.xnames=arrset(self.xnames,iv,'A2')
        [iv,ix_,_] = s2mpj_ii('A4',ix_)
        self.xnames=arrset(self.xnames,iv,'A4')
        [iv,ix_,_] = s2mpj_ii('B',ix_)
        self.xnames=arrset(self.xnames,iv,'B')
        [iv,ix_,_] = s2mpj_ii('C',ix_)
        self.xnames=arrset(self.xnames,iv,'C')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['1']),int(v_['M'])+1):
            v_['XSQR'] = v_['X'+str(I)]*v_['X'+str(I)]
            v_['XQUART'] = v_['XSQR']*v_['XSQR']
            [ig,ig_,_] = s2mpj_ii('O'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['A2']])
            valA = np.append(valA,float(v_['XSQR']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['A4']])
            valA = np.append(valA,float(v_['XQUART']))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        selfnob      = ngrp
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        for I in range(int(v_['1']),int(v_['M'])+1):
            self.gconst = arrset(self.gconst,ig_['O'+str(I)],float(v_['Y'+str(I)]))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        self.xlower[ix_['A2']] = -float('Inf')
        self.xupper[ix_['A2']] = +float('Inf')
        self.xlower[ix_['A4']] = -float('Inf')
        self.xupper[ix_['A4']] = +float('Inf')
        self.xlower[ix_['B']] = 0.00001
        self.xlower[ix_['C']] = 0.00001
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(1.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eQUOT', iet_)
        elftv = loaset(elftv,it,0,'B')
        elftv = loaset(elftv,it,1,'C')
        elftp = []
        elftp = loaset(elftp,it,0,'XSQR')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['1']),int(v_['M'])+1):
            v_['XSQR'] = v_['X'+str(I)]*v_['X'+str(I)]
            ename = 'E'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eQUOT')
            ielftype = arrset(ielftype,ie,iet_["eQUOT"])
            vname = 'B'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(1.0))
            posev = np.where(elftv[ielftype[ie]]=='B')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'C'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(1.0))
            posev = np.where(elftv[ielftype[ie]]=='C')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='XSQR')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['XSQR']))
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
        for I in range(int(v_['1']),int(v_['M'])+1):
            ig = ig_['O'+str(I)]
            self.grftype = arrset(self.grftype,ig,'gL2')
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        self.objlower = 0.0
#    Solution
# LO SOLTN                6.83513867
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass   = "C-CSBR2-RN-4-0"
        self.objderlvl = 2

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eQUOT(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        DENOM = 1.0/(EV_[1,0]+self.elpar[iel_][0])
        f_   = EV_[0,0]*DENOM
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = DENOM
            g_[1] = -EV_[0,0]*DENOM*DENOM
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = -DENOM*DENOM
                H_[1,0] = H_[0,1]
                H_[1,1] = 2.0*EV_[0,0]*DENOM**3
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

