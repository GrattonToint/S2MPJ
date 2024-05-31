from s2xlib import *
class  READING2(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : READING2
#    *********
# 
#    A linear optimal control problem from Nancy Nichols
#    with a given initial condition.
#    This problem arises in tide modelling.
# 
#    Source:
#    S. Lyle and N.K. Nichols,
#    "Numerical Methods for Optimal Control Problems with State Constraints",
#    Numerical Analysis Report 8/91, Dept of Mathematics, 
#    University of Reading, UK.
# 
#    SIF input: Nick Gould, July 1991.
# 
#    classification = "LLR2-MN-V-V"
# 
#    Number of discretized points in [0,1]
# 
#           Alternative values for the SIF file parameters:
# IE N                   2              $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'READING2'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'READING2'
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
# IE N                   50             $-PARAMETER
# IE N                   100            $-PARAMETER     original value
# IE N                   500            $-PARAMETER
# IE N                   1000           $-PARAMETER
# IE N                   2000           $-PARAMETER
# IE N                   5000           $-PARAMETER
        v_['PI'] = 3.1415926535
        v_['2PI'] = 2.0*v_['PI']
        v_['PI**2'] = v_['PI']*v_['PI']
        v_['8PI**2'] = 8.0*v_['PI**2']
        v_['1/8PI**2'] = 1.0/v_['8PI**2']
        v_['A'] = 0.07716
        v_['1/A'] = 1.0/v_['A']
        v_['1/2A'] = 2.0*v_['1/A']
        v_['2A'] = 2.0*v_['A']
        v_['-2A'] = -1.0*v_['2A']
        v_['-1/2A'] = 1.0/v_['-2A']
        v_['N-1'] = -1+v_['N']
        v_['RN'] = float(v_['N'])
        v_['H'] = 1.0/v_['RN']
        v_['2/H'] = 2.0*v_['RN']
        v_['H/2'] = 0.5*v_['H']
        v_['-H/2'] = -1.0*v_['H/2']
        v_['1/H'] = 1.0*v_['RN']
        v_['-1/H'] = -1.0*v_['RN']
        v_['H/8PI**2'] = v_['1/8PI**2']*v_['H']
        v_['0'] = 0
        v_['1'] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['0']),int(v_['N'])+1):
            [iv,ix_,_] = s2x_ii('X1u'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X1u'+str(I))
            [iv,ix_,_] = s2x_ii('X2u'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X2u'+str(I))
            [iv,ix_,_] = s2x_ii('U'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'U'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['RI'] = float(I)
            v_['TI'] = v_['RI']*v_['H']
            v_['2PITI'] = v_['2PI']*v_['TI']
            v_['CTI'] = np.cos(v_['2PITI'])
            v_['-CCTI'] = v_['CTI']*v_['-H/2']
            v_['I-1'] = -1+I
            v_['RI-1'] = float(v_['I-1'])
            v_['TI-1'] = v_['RI-1']*v_['H']
            v_['2PITI-1'] = v_['2PI']*v_['TI-1']
            v_['CTI-1'] = np.cos(v_['2PITI-1'])
            v_['-CCTI-1'] = v_['CTI-1']*v_['-H/2']
            [ig,ig_,_] = s2x_ii('COST',ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X1u'+str(I)]
            pbm.A[ig,iv] = float(v_['-CCTI'])+pbm.A[ig,iv]
            iv = ix_['X1u'+str(int(v_['I-1']))]
            pbm.A[ig,iv] = float(v_['-CCTI-1'])+pbm.A[ig,iv]
            iv = ix_['U'+str(I)]
            pbm.A[ig,iv] = float(v_['H/8PI**2'])+pbm.A[ig,iv]
            iv = ix_['U'+str(int(v_['I-1']))]
            pbm.A[ig,iv] = float(v_['H/8PI**2'])+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['I-1'] = -1+I
            [ig,ig_,_] = s2x_ii('C1u'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'C1u'+str(I))
            iv = ix_['X1u'+str(I)]
            pbm.A[ig,iv] = float(v_['1/H'])+pbm.A[ig,iv]
            iv = ix_['X1u'+str(int(v_['I-1']))]
            pbm.A[ig,iv] = float(v_['-1/H'])+pbm.A[ig,iv]
            iv = ix_['X2u'+str(I)]
            pbm.A[ig,iv] = float(-0.5)+pbm.A[ig,iv]
            iv = ix_['X2u'+str(int(v_['I-1']))]
            pbm.A[ig,iv] = float(-0.5)+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('C2u'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'C2u'+str(I))
            iv = ix_['X2u'+str(I)]
            pbm.A[ig,iv] = float(v_['1/H'])+pbm.A[ig,iv]
            iv = ix_['X2u'+str(int(v_['I-1']))]
            pbm.A[ig,iv] = float(v_['-1/H'])+pbm.A[ig,iv]
            iv = ix_['U'+str(I)]
            pbm.A[ig,iv] = float(-0.5)+pbm.A[ig,iv]
            iv = ix_['U'+str(int(v_['I-1']))]
            pbm.A[ig,iv] = float(-0.5)+pbm.A[ig,iv]
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
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('inf'))
        pb.xupper = np.full((pb.n,1),float('inf'))
        pb.xlower[ix_['X1u'+str(int(v_['0']))]] = 0.0
        pb.xupper[ix_['X1u'+str(int(v_['0']))]] = 0.0
        pb.xlower[ix_['X2u'+str(int(v_['0']))]] = 0.0
        pb.xupper[ix_['X2u'+str(int(v_['0']))]] = 0.0
        for I in range(int(v_['1']),int(v_['N'])+1):
            pb.xlower[ix_['X1u'+str(I)]] = -float('Inf')
            pb.xupper[ix_['X1u'+str(I)]] = +float('Inf')
            pb.xlower[ix_['X2u'+str(I)]] = -0.125
            pb.xupper[ix_['X2u'+str(I)]] = 0.125
        for I in range(int(v_['0']),int(v_['N'])+1):
            pb.xlower[ix_['U'+str(I)]] = -1.0
            pb.xupper[ix_['U'+str(I)]] = 1.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.lincons   = np.arange(len(pbm.congrps))
        pb.pbclass = "LLR2-MN-V-V"
        pb.x0          = np.zeros((pb.n,1))
        self.pb = pb; self.pbm = pbm

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

