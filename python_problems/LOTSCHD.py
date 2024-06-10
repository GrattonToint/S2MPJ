from s2mpjlib import *
class  LOTSCHD(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    A simple quadratic program inspired by the economic lot scheduling
#    problem.
# 
#    Source:
#    an exercize for L. Watson course on LANCELOT in the Spring 1993.
# 
#    SIF input: T. Kuan, Virginia Tech., Spring 1993.
# 
#    classification = "QLR2-AN-12-7"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'LOTSCHD'

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
        v_['X1'] = 1.502
        v_['X2'] = 1.126
        v_['X3'] = 0.815
        v_['X4'] = 1.268
        v_['X5'] = 1.502
        v_['X6'] = 0.740
        v_['A1'] = 1.8
        v_['A2'] = 3.2
        v_['A3'] = 6.1
        v_['A4'] = 3.2
        v_['A5'] = 1.8
        v_['A6'] = 7.4
        v_['C1'] = 11.0
        v_['C2'] = 3.0
        v_['C3'] = 20.0
        v_['C4'] = 17.0
        v_['C5'] = 9.0
        v_['C6'] = 20.0
        v_['C7'] = 126.1
        v_['1'] = 1
        v_['2'] = 2
        v_['4'] = 4
        v_['5'] = 5
        v_['6'] = 6
        v_['7'] = 7
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['6'])+1):
            [iv,ix_,_] = s2mpj_ii('T'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'T'+str(I))
            [iv,ix_,_] = s2mpj_ii('U'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'U'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['6'])+1):
            [ig,ig_,_] = s2mpj_ii('OBJ'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['T'+str(I)]
            pbm.A[ig,iv] = float(v_['X'+str(I)])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('CONS7',ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'CONS7')
            iv = ix_['T'+str(I)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['U'+str(I)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('CONS'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'CONS'+str(I))
            iv = ix_['T'+str(I)]
            pbm.A[ig,iv] = float(v_['A'+str(I)])+pbm.A[ig,iv]
            iv = ix_['U'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        for I in range(int(v_['2']),int(v_['4'])+1):
            [ig,ig_,_] = s2mpj_ii('CONS'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'CONS'+str(I))
            iv = ix_['T'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            iv = ix_['U'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('CONS2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'CONS2')
        iv = ix_['T3']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['U3']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['2'])+1):
            [ig,ig_,_] = s2mpj_ii('CONS3',ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'CONS3')
            iv = ix_['T'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            iv = ix_['U'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        for I in range(int(v_['4']),int(v_['6'])+1):
            [ig,ig_,_] = s2mpj_ii('CONS3',ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'CONS3')
            iv = ix_['T'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            iv = ix_['U'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('CONS4',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'CONS4')
        iv = ix_['T1']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['U1']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        for I in range(int(v_['5']),int(v_['6'])+1):
            [ig,ig_,_] = s2mpj_ii('CONS4',ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'CONS4')
            iv = ix_['T'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            iv = ix_['U'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('CONS5',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'CONS5')
        iv = ix_['T1']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['U1']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['5'])+1):
            [ig,ig_,_] = s2mpj_ii('CONS6',ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'CONS6')
            iv = ix_['T'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            iv = ix_['U'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
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
        for I in range(int(v_['1']),int(v_['7'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['CONS'+str(I)],float(v_['C'+str(I)]))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gSQUARE',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['6'])+1):
            ig = ig_['OBJ'+str(I)]
            pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
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
        pb.pbclass = "QLR2-AN-12-7"
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

