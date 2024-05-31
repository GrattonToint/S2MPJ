from s2xlib import *
class  DEGTRID(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DEGTRID
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

    name = 'DEGTRID'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'DEGTRID'
        pbm.name = self.name
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
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['0']),int(v_['N'])+1):
            [iv,ix_,_] = s2x_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2x_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X'+str(int(v_['0']))]
        pbm.A[ig,iv] = float(-0.5)+pbm.A[ig,iv]
        iv = ix_['X'+str(int(v_['1']))]
        pbm.A[ig,iv] = float(-1.5)+pbm.A[ig,iv]
        for I in range(int(v_['2']),int(v_['N-1'])+1):
            [ig,ig_,_] = s2x_ii('OBJ',ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = float(-2.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X'+str(int(v_['N']))]
        pbm.A[ig,iv] = float(-1.5)+pbm.A[ig,iv]
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('inf'))
        pb.xupper = np.full((pb.n,1),float('inf'))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(2.0))
        #%%%%%%%%%%%%%%%%%%%% QUADRATIC %%%%%%%%%%%%%%%%%%%
        pbm.H = lil_matrix((pb.n, pb.n))
        ix1 = ix_['X'+str(int(v_['0']))]
        ix2 = ix_['X'+str(int(v_['0']))]
        pbm.H[ix1,ix2] = float(1.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['I-1'] = -1+I
            ix1 = ix_['X'+str(I)]
            ix2 = ix_['X'+str(I)]
            pbm.H[ix1,ix2] = float(1.0)+pbm.H[ix1,ix2]
            pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
            ix1 = ix_['X'+str(I)]
            ix2 = ix_['X'+str(int(v_['I-1']))]
            pbm.H[ix1,ix2] = float(0.5)+pbm.H[ix1,ix2]
            pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "QBR2-AN-V-0"
        pbm.H = pbm.H.tocsr()
        self.pb = pb; self.pbm = pbm

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

