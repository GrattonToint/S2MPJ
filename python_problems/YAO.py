from s2mpjlib import *
class  YAO(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
#    A linear least-sqaure problem with k-convex constraints
#       min (1/2) || f(t) - x ||^2
#    subject to the constraints
#       _ 
#       \/_k  x  >=  0,
#    where  f(t) and  x  are vectors in (n+k)-dimensional space.
# 
#    We choose f(t) = sin(t), x(1) >= 0.08 and fix x(n+i) = 0
# 
#    SIF input: Aixiang Yao, Virginia Tech., May 1995
#               modifications by Nick Gould
# 
#    classification = "QLR2-AN-V-V"
# 
#   Number of discretization points
# 
#           Alternative values for the SIF file parameters:
# IE P                   20             $-PARAMETER
# IE P                   200            $-PARAMETER
# IE P                   2000           $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'YAO'

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
        if nargin<1:
            v_['P'] = int(20);  #  SIF file default value
        else:
            v_['P'] = int(args[0])
# IE k                   2              $-PARAMETER
        if nargin<2:
            v_['k'] = int(2);  #  SIF file default value
        else:
            v_['k'] = int(args[1])
# IE k                   3              $-PARAMETER
# IE k                   4              $-PARAMETER
        v_['1'] = 1
        v_['2'] = 2
        v_['P+1'] = v_['P']+v_['1']
        v_['P+k'] = v_['P']+v_['k']
        v_['RP'] = float(v_['P+k'])
        v_['OVP'] = 1.0/v_['RP']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for i in range(int(v_['1']),int(v_['P+k'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(i),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(i))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for i in range(int(v_['1']),int(v_['P+k'])+1):
            [ig,ig_,_] = s2mpj_ii('S'+str(i),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(i)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            pbm.gscale = arrset(pbm.gscale,ig,float(2.0))
        for i in range(int(v_['1']),int(v_['P'])+1):
            v_['i+1'] = 1+i
            [ig,ig_,_] = s2mpj_ii('B'+str(i),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'B'+str(i))
            iv = ix_['X'+str(i)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['i+1']))]
            pbm.A[ig,iv] = float(-2.0)+pbm.A[ig,iv]
            v_['i+2'] = 2+i
            iv = ix_['X'+str(int(v_['i+2']))]
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
        for i in range(int(v_['1']),int(v_['P+k'])+1):
            v_['Ri'] = float(i)
            v_['iOVP'] = v_['Ri']*v_['OVP']
            v_['SINI'] = np.sin(v_['iOVP'])
            pbm.gconst = arrset(pbm.gconst,ig_['S'+str(i)],float(v_['SINI']))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower = np.zeros((pb.n,1))
        pb.xlower[ix_['X'+str(int(v_['1']))]] = 0.08
        for i in range(int(v_['P+1']),int(v_['P+k'])+1):
            pb.xlower[ix_['X'+str(i)]] = 0.0
            pb.xupper[ix_['X'+str(i)]] = 0.0
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gSQ',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for i in range(int(v_['1']),int(v_['P+k'])+1):
            ig = ig_['S'+str(i)]
            pbm.grftype = arrset(pbm.grftype,ig,'gSQ')
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# XL SOLUTION             2.39883D+00   $ (p=20)
# XL SOLUTION             2.01517D+01   $ (p=200)
# XL SOLUTION             1.97705D+02   $ (p=2000)
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
        pb.pbclass = "QLR2-AN-V-V"
        pb.x0          = np.zeros((pb.n,1))
        self.pb = pb; self.pbm = pbm
# ********************
#  SET UP THE GROUPS *
#  ROUTINE           *
# ********************

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gSQ(pbm,nargout,*args):

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

