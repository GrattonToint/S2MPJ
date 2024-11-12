from s2mpjlib import *
class  ALSOTAME(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : ALSOTAME
#    *********
#    Another simple constrained  problem
# 
#    Source:
#    A.R. Conn, N. Gould and Ph.L. Toint,
#    "The LANCELOT User's Manual",
#    Dept of Maths, FUNDP, 1991.
# 
#    SIF input:  Ph. Toint, Jan 1991.
# 
#    classification = "C-COOR2-AN-2-1"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 9 XI 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'ALSOTAME'

    def __init__(self, *args): 
        import numpy as np
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2mpj_ii('x',ix_)
        self.xnames=arrset(self.xnames,iv,'x')
        [iv,ix_,_] = s2mpj_ii('y',ix_)
        self.xnames=arrset(self.xnames,iv,'y')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.A       = lil_matrix((1000000,1000000))
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames      = np.array([])
        self.cnames = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('Object',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['x']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['y']
        self.A[ig,iv] = float(-2.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('Constr',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'Constr')
        iv = ix_['x']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['y']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        legrps = np.where(gtype=='<=')[0]
        eqgrps = np.where(gtype=='==')[0]
        gegrps = np.where(gtype=='>=')[0]
        self.nle = len(legrps)
        self.neq = len(eqgrps)
        self.nge = len(gegrps)
        self.m   = self.nle+self.neq+self.nge
        self.congrps = np.concatenate((legrps,eqgrps,gegrps))
        self.cnames= cnames[self.congrps]
        self.nob = ngrp-self.m
        self.objgrps = np.where(gtype=='<>')[0]
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        self.gconst = arrset(self.gconst,ig_['Constr'],float(1.0))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        self.xlower[ix_['x']] = -2.0
        self.xupper[ix_['x']] = 2.0
        self.xlower[ix_['y']] = -1.5
        self.xupper[ix_['y']] = 1.5
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gEXPN',igt_)
        [it,igt_,_] = s2mpj_ii('gSINE',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        ig = ig_['Object']
        self.grftype = arrset(self.grftype,ig,'gEXPN')
        ig = ig_['Constr']
        self.grftype = arrset(self.grftype,ig,'gSINE')
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        self.A.resize(ngrp,self.n)
        self.A     = self.A.tocsr()
        sA1,sA2    = self.A.shape
        self.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons   = np.arange(len(self.congrps))
        self.pbclass = "C-COOR2-AN-2-1"
        self.x0        = np.zeros((self.n,1))
        self.objderlvl = 2
        self.conderlvl = [2]


    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gEXPN(self,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        EXPA = np.exp(GVAR_)
        f_= EXPA
        if nargout>1:
            g_ = EXPA
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = EXPA
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def gSINE(self,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        SINA = np.sin(GVAR_)
        f_= SINA
        if nargout>1:
            g_ = np.cos(GVAR_)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = -SINA
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

