from s2mpjlib import *
class  DEGTRID2(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DEGTRID2
#    *********
# 
#    A degenerate bound constrained convex quadratic program
#    with a tri-diagonal Hessian
# 
#    SIF input: Nick Gould, August 2011
# 
#    classification = "QBR2-AN-V-0"
# 
#    The number of variables - 1
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'DEGTRID2'

    def __init__(self, *args): 
        import numpy as np
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(10);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
#           Alternative values for the SIF file parameters:
# IE N                   50             $-PARAMETER
# IE N                   100            $-PARAMETER
# IE N                   1000           $-PARAMETER
# IE N                   10000          $-PARAMETER
# IE N                   100000         $-PARAMETER
        v_['0'] = 0
        v_['1'] = 1
        v_['2'] = 2
        v_['N-1'] = -1+v_['N']
        v_['N+1'] = 1+v_['N']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['0']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.A       = lil_matrix((1000000,1000000))
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames      = np.array([])
        self.cnames = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X'+str(int(v_['0']))]
        self.A[ig,iv] = float(-0.5)+self.A[ig,iv]
        iv = ix_['X'+str(int(v_['1']))]
        self.A[ig,iv] = float(-1.5)+self.A[ig,iv]
        for I in range(int(v_['2']),int(v_['N-1'])+1):
            [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(I)]
            self.A[ig,iv] = float(-2.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X'+str(int(v_['N']))]
        self.A[ig,iv] = float(-1.5)+self.A[ig,iv]
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),1.0)
        self.xupper = np.full((self.n,1),+float('inf'))
        self.xlower[ix_['X0']] = 0.0
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(2.0))
        #%%%%%%%%%%%%%%%%%%%% QUADRATIC %%%%%%%%%%%%%%%%%%%
        self.H = lil_matrix((self.n, self.n))
        ix1 = ix_['X'+str(int(v_['0']))]
        ix2 = ix_['X'+str(int(v_['0']))]
        self.H[ix1,ix2] = float(1.0)+self.H[ix1,ix2]
        self.H[ix2,ix1] = self.H[ix1,ix2]
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['I-1'] = -1+I
            ix1 = ix_['X'+str(I)]
            ix2 = ix_['X'+str(I)]
            self.H[ix1,ix2] = float(1.0)+self.H[ix1,ix2]
            self.H[ix2,ix1] = self.H[ix1,ix2]
            ix1 = ix_['X'+str(I)]
            ix2 = ix_['X'+str(int(v_['I-1']))]
            self.H[ix1,ix2] = float(0.5)+self.H[ix1,ix2]
            self.H[ix2,ix1] = self.H[ix1,ix2]
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        self.A.resize(ngrp,self.n)
        self.A     = self.A.tocsr()
        sA1,sA2    = self.A.shape
        self.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass = "QBR2-AN-V-0"
        self.H = self.H.tocsr()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

