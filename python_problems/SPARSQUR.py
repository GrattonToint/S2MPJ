from s2mpjlib import *
class  SPARSQUR(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : SPARSQUR
#    *********
# 
#    A sparse quartic problem
# 
#    SIF input: Nick Gould, November 1995
# 
#    classification = "C-COUR2-AN-V-0"
# 
#    The number of variables 
# 
#           Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER
# IE N                   50             $-PARAMETER
# IE N                   100            $-PARAMETER
# IE N                   1000           $-PARAMETER     original value
# IE N                   5000           $-PARAMETER
# IE N                   10000          $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'SPARSQUR'

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
        v_['1'] = 1
        v_['2'] = 2
        v_['4'] = 4
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
        for I in range(int(v_['1']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('OBJ'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        selfnob      = ngrp
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(0.5))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQUR', iet_)
        elftv = loaset(elftv,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['1']),int(v_['N'])+1):
            ename = 'S'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSQUR')
            ielftype = arrset(ielftype,ie,iet_["eSQUR"])
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.5))
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gSQR',igt_)
        [it,igt_,_] = s2mpj_ii('gSQR',igt_)
        grftp = []
        grftp = loaset(grftp,it,0,'P')
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        self.grpar   = []
        for I in range(int(v_['1']),int(v_['N'])+1):
            ig = ig_['OBJ'+str(I)]
            self.grftype = arrset(self.grftype,ig,'gSQR')
            v_['RI'] = float(I)
            posgp = np.where(grftp[igt_[self.grftype[ig]]]=='P')[0]
            self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['RI']))
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['S'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,float(1.0))
            v_['J'] = 2*I
            v_['J'] = -1+v_['J']
            v_['K'] = int(np.fix(v_['J']/v_['N']))
            v_['K'] = v_['K']*v_['N']
            v_['J'] = v_['J']-v_['K']
            v_['J'] = 1+v_['J']
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['S'+str(int(v_['J']))])
            self.grelw = loaset(self.grelw,ig,posel,float(1.0))
            v_['J'] = 3*I
            v_['J'] = -1+v_['J']
            v_['K'] = int(np.fix(v_['J']/v_['N']))
            v_['K'] = v_['K']*v_['N']
            v_['J'] = v_['J']-v_['K']
            v_['J'] = 1+v_['J']
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['S'+str(int(v_['J']))])
            self.grelw = loaset(self.grelw,ig,posel,float(1.0))
            v_['J'] = 5*I
            v_['J'] = -1+v_['J']
            v_['K'] = int(np.fix(v_['J']/v_['N']))
            v_['K'] = v_['K']*v_['N']
            v_['J'] = v_['J']-v_['K']
            v_['J'] = 1+v_['J']
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['S'+str(int(v_['J']))])
            self.grelw = loaset(self.grelw,ig,posel,float(1.0))
            v_['J'] = 7*I
            v_['J'] = -1+v_['J']
            v_['K'] = int(np.fix(v_['J']/v_['N']))
            v_['K'] = v_['K']*v_['N']
            v_['J'] = v_['J']-v_['K']
            v_['J'] = 1+v_['J']
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['S'+str(int(v_['J']))])
            self.grelw = loaset(self.grelw,ig,posel,float(1.0))
            v_['J'] = 11*I
            v_['J'] = -1+v_['J']
            v_['K'] = int(np.fix(v_['J']/v_['N']))
            v_['K'] = v_['K']*v_['N']
            v_['J'] = v_['J']-v_['K']
            v_['J'] = 1+v_['J']
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['S'+str(int(v_['J']))])
            self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLUTION            0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass   = "C-COUR2-AN-V-0"
        self.objderlvl = 2

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSQUR(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = 0.5e0*EV_[0,0]*EV_[0,0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[0,0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 1.0e0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gSQR(self,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= 0.5e0*self.grpar[igr_][0]*GVAR_*GVAR_
        if nargout>1:
            g_ = self.grpar[igr_][0]*GVAR_
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = self.grpar[igr_][0]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

