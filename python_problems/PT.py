from s2mpjlib import *
class  PT(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : PT
#    *********
# 
#    A nonlinear programming formulation of a discretization of
#    a nonlinear minimax problem.
# 
#    The problem is
# 
#        min  max phi(x,w), for all w in the interval I.
#         x    w
# 
#    I is discretized, and the problem solved over the
#    discrete points.
# 
#    Nonlinear programming formulation
#        min   u     s.t.  u - phi >= 0
#        x,u
# 
#    Specific problem: I = [0,1]
#    phi(x,w) = (2w^2-1)x + w(1-w)(1-x)
# 
#    Source: E. Polak and A. L. Tits,
#    "A recursive quadratic programming algorithm for semi-infinite
#     optimization problems",
#    Appl. Math. Optim. 8, 1982, pp 325-349.
# 
#    SIF input: Nick Gould, February, 1994.
# 
#    classification = "C-CLLR2-AN-2-V"
# 
#    Discretization
# 
# IE M                   2
# IE M                   100
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'PT'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['M'] = 500
        v_['LOWER'] = 0.0
        v_['UPPER'] = 1.0
        v_['ONE'] = 1.0
        v_['0'] = 0
        v_['DIFF'] = v_['UPPER']-v_['LOWER']
        v_['RM'] = float(v_['M'])
        v_['H'] = v_['DIFF']/v_['RM']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [iv,ix_,_] = s2mpj_ii('U',ix_)
        self.xnames=arrset(self.xnames,iv,'U')
        [iv,ix_,_] = s2mpj_ii('X',ix_)
        self.xnames=arrset(self.xnames,iv,'X')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['U']])
        valA = np.append(valA,float(1.0))
        for I in range(int(v_['0']),int(v_['M'])+1):
            v_['RI'] = float(I)
            v_['W'] = v_['RI']*v_['H']
            v_['W'] = v_['W']+v_['LOWER']
            v_['1-W'] = 1.0-v_['W']
            v_['W(1-W)'] = v_['W']*v_['1-W']
            v_['W**2'] = v_['W']*v_['W']
            v_['2W**2'] = 2.0*v_['W**2']
            v_['2W**2-1'] = v_['2W**2']-v_['ONE']
            v_['XCOEFF'] = v_['W(1-W)']-v_['2W**2-1']
            [ig,ig_,_] = s2mpj_ii('LO'+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'LO'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['U']])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X']])
            valA = np.append(valA,float(v_['XCOEFF']))
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
        self.cnames = cnames[self.congrps]
        self.nob = ngrp-self.m
        self.objgrps = np.where(gtype=='<>')[0]
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        for I in range(int(v_['0']),int(v_['M'])+1):
            v_['RI'] = float(I)
            v_['W'] = v_['RI']*v_['H']
            v_['W'] = v_['W']+v_['LOWER']
            v_['1-W'] = 1.0-v_['W']
            v_['W-W**2'] = v_['W']*v_['1-W']
            self.gconst = arrset(self.gconst,ig_['LO'+str(I)],float(v_['W-W**2']))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle+self.neq,self.m)] = np.zeros((self.nge,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons   = np.arange(len(self.congrps))
        self.pbclass   = "C-CLLR2-AN-2-V"
        self.x0        = np.zeros((self.n,1))
        self.objderlvl = 2
        self.conderlvl = [2]


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

