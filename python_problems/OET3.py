from s2mpjlib import *
class  OET3(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : OET3
#    *********
# 
#    A nonlinear programming formulation of a discretization of
#    a nonlinear Chebychev problem.
# 
#    The problem is
# 
#        min  max || phi(x,w) ||, for all w in the interval I.
#         x    w
# 
#    I is discretized, and the problem solved over the
#    discrete points.
# 
#    Nonlinear programming formulation
#        min   u     s.t.  u - phi >= 0, u + phi >= 0
#        x,u
# 
#    Specific problem: I = [0,1]
#    phi(x,w) = sin(w) - (x1 + x2 w + x3 w^2 )
# 
#    Source: K. Oettershagen
#    "Ein superlinear konvergenter algorithmus zur losung 
#     semi-infiniter optimierungsproblem",
#     Ph.D thesis, Bonn University, 1982
# 
#    SIF input: Nick Gould, February, 1994.
# 
#    classification = "LLR2-AN-4-V"
# 
#    Discretization
# 
# IE M                   2
# IE M                   100
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'OET3'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['M'] = 500
        v_['LOWER'] = 0.0
        v_['UPPER'] = 1.0
        v_['0'] = 0
        v_['DIFF'] = v_['UPPER']-v_['LOWER']
        v_['RM'] = float(v_['M'])
        v_['H'] = v_['DIFF']/v_['RM']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2mpj_ii('U',ix_)
        pb.xnames=arrset(pb.xnames,iv,'U')
        [iv,ix_,_] = s2mpj_ii('X1',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X1')
        [iv,ix_,_] = s2mpj_ii('X2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X2')
        [iv,ix_,_] = s2mpj_ii('X3',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X3')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['U']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        for I in range(int(v_['0']),int(v_['M'])+1):
            v_['RI'] = float(I)
            v_['W'] = v_['RI']*v_['H']
            v_['W'] = v_['W']+v_['LOWER']
            v_['-W'] = -1.0*v_['W']
            v_['W**2'] = v_['W']*v_['W']
            v_['-W**2'] = -1.0*v_['W**2']
            [ig,ig_,_] = s2mpj_ii('LO'+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'LO'+str(I))
            iv = ix_['U']
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['X1']
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            iv = ix_['X2']
            pbm.A[ig,iv] = float(v_['-W'])+pbm.A[ig,iv]
            iv = ix_['X3']
            pbm.A[ig,iv] = float(v_['-W**2'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('UP'+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'UP'+str(I))
            iv = ix_['U']
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['X1']
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['X2']
            pbm.A[ig,iv] = float(v_['W'])+pbm.A[ig,iv]
            iv = ix_['X3']
            pbm.A[ig,iv] = float(v_['W**2'])+pbm.A[ig,iv]
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        legrps = find(gtype,lambda x:x=='<=')
        eqgrps = find(gtype,lambda x:x=='==')
        gegrps = find(gtype,lambda x:x=='>=')
        pb.nle = len(legrps)
        pb.neq = len(eqgrps)
        pb.nge = len(gegrps)
        pb.m   = pb.nle+pb.neq+pb.nge
        pbm.congrps = find(gtype,lambda x:(x=='<=' or x=='==' or x=='>='))
        pb.cnames= cnames[pbm.congrps]
        pb.nob = ngrp-pb.m
        pbm.objgrps = find(gtype,lambda x:x=='<>')
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        for I in range(int(v_['0']),int(v_['M'])+1):
            v_['RI'] = float(I)
            v_['W'] = v_['RI']*v_['H']
            v_['W'] = v_['W']+v_['LOWER']
            v_['SINW'] = np.sin(v_['W'])
            v_['-SINW'] = -1.0*v_['SINW']
            pbm.gconst = arrset(pbm.gconst,ig_['LO'+str(I)],float(v_['-SINW']))
            pbm.gconst = arrset(pbm.gconst,ig_['UP'+str(I)],float(v_['SINW']))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower = np.zeros((pb.n,1))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle+pb.neq,pb.m)] = np.zeros((pb.nge,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.lincons   = np.arange(len(pbm.congrps))
        pb.pbclass = "LLR2-AN-4-V"
        pb.x0          = np.zeros((pb.n,1))
        self.pb = pb; self.pbm = pbm

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

