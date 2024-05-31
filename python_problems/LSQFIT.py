from s2xlib import *
class  LSQFIT(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : LSQFIT
#    *********
#    An elementary constrained linear least-squares fit
# 
#    Source:
#    A.R. Conn, N. Gould and Ph.L. Toint,
#    "The LANCELOT User's Manual",
#    Dept of Maths, FUNDP, 1991.
# 
#    SIF input: Ph. Toint, Jan 1991.
# 
#    classification = "SLR2-AN-2-1"
# 
#    Data points
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'LSQFIT'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'LSQFIT'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['X1'] = 0.1
        v_['X2'] = 0.3
        v_['X3'] = 0.5
        v_['X4'] = 0.7
        v_['X5'] = 0.9
        v_['Y1'] = 0.25
        v_['Y2'] = 0.3
        v_['Y3'] = 0.625
        v_['Y4'] = 0.701
        v_['Y5'] = 1.0
        v_['C'] = 0.85
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2x_ii('a',ix_)
        pb.xnames=arrset(pb.xnames,iv,'a')
        [iv,ix_,_] = s2x_ii('b',ix_)
        pb.xnames=arrset(pb.xnames,iv,'b')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2x_ii('Obj1',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['a']
        pbm.A[ig,iv] = float(v_['X1'])+pbm.A[ig,iv]
        iv = ix_['b']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(2.0))
        [ig,ig_,_] = s2x_ii('Obj2',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['a']
        pbm.A[ig,iv] = float(v_['X2'])+pbm.A[ig,iv]
        iv = ix_['b']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(2.0))
        [ig,ig_,_] = s2x_ii('Obj3',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['a']
        pbm.A[ig,iv] = float(v_['X3'])+pbm.A[ig,iv]
        iv = ix_['b']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(2.0))
        [ig,ig_,_] = s2x_ii('Obj4',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['a']
        pbm.A[ig,iv] = float(v_['X4'])+pbm.A[ig,iv]
        iv = ix_['b']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(2.0))
        [ig,ig_,_] = s2x_ii('Obj5',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['a']
        pbm.A[ig,iv] = float(v_['X5'])+pbm.A[ig,iv]
        iv = ix_['b']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(2.0))
        [ig,ig_,_] = s2x_ii('Cons',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'Cons')
        iv = ix_['a']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['b']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
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
        pbm.gconst = arrset(pbm.gconst,ig_['Obj1'],float(v_['Y1']))
        pbm.gconst = arrset(pbm.gconst,ig_['Obj2'],float(v_['Y2']))
        pbm.gconst = arrset(pbm.gconst,ig_['Obj3'],float(v_['Y3']))
        pbm.gconst = arrset(pbm.gconst,ig_['Obj4'],float(v_['Y4']))
        pbm.gconst = arrset(pbm.gconst,ig_['Obj5'],float(v_['Y5']))
        pbm.gconst = arrset(pbm.gconst,ig_['Cons'],float(v_['C']))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('inf'))
        pb.xupper = np.full((pb.n,1),float('inf'))
        pb.xlower[ix_['b']] = -float('Inf')
        pb.xupper[ix_['b']] = +float('Inf')
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2x_ii('gSQUARE',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        ig = ig_['Obj1']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        ig = ig_['Obj2']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        ig = ig_['Obj3']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        ig = ig_['Obj4']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        ig = ig_['Obj5']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.cupper[np.arange(pb.nle)] = np.zeros((pb.nle,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.lincons   = np.arange(len(pbm.congrps))
        pb.pbclass = "SLR2-AN-2-1"
        pb.x0          = np.zeros((pb.n,1))
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gSQUARE(pbm,nargout,*args):

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

