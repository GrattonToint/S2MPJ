from s2xlib import *
class  REPEAT(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : REPEAT
#    *********
# 
#    This problem is to find the nearest feasible point to 2n+1 inconsistent
#    linear equations subject to bounds
# 
#    Source: blue-cheese delerium
# 
#    SIF input: Nick Gould, December 2020.
# 
#    classification = "NLR2-AN-V-V"
# 
#    N is the number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'REPEAT'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'REPEAT'
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
        v_['N-1'] = -1+v_['N']
        v_['1'] = 1
        v_['2'] = 2
        v_['100'] = 100
        v_['N/2'] = int(np.fix(v_['N']/v_['2']))
        v_['N/100'] = int(np.fix(v_['N']/v_['100']))
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2x_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            v_['I+1'] = 1+I
            [ig,ig_,_] = s2x_ii('C'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'C'+str(I))
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = 1.0+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['I+1']))]
            pbm.A[ig,iv] = 1.0+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('C'+str(int(v_['N'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C'+str(int(v_['N'])))
        iv = ix_['X'+str(int(v_['N']))]
        pbm.A[ig,iv] = 1.0+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            v_['I+1'] = 1+I
            [ig,ig_,_] = s2x_ii('R'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'R'+str(I))
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = 1.0+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['I+1']))]
            pbm.A[ig,iv] = 1.0+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('R'+str(int(v_['N'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'R'+str(int(v_['N'])))
        iv = ix_['X'+str(int(v_['N']))]
        pbm.A[ig,iv] = 1.0+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['N'])+1,int(v_['N/100'])):
            v_['RI'] = float(I)
            [ig,ig_,_] = s2x_ii('E',ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'E')
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = v_['RI']+pbm.A[ig,iv]
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
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['C'+str(I)],2.0)
        pbm.gconst = arrset(pbm.gconst,ig_['C'+str(int(v_['N']))],1.0)
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['R'+str(I)],4.0)
        pbm.gconst = arrset(pbm.gconst,ig_['R'+str(int(v_['N']))],3.0)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower = np.full((pb.n,1),-1.0)
        pb.xupper[ix_['X'+str(int(v_['1']))]] = 0.0
        pb.xlower[ix_['X'+str(int(v_['2']))]] = 3.0
        pb.xupper[ix_['X'+str(int(v_['N/2']))]] = 0.0
        pb.xlower[ix_['X'+str(int(v_['N-1']))]] = 3.0
        pb.xupper[ix_['X'+str(int(v_['N']))]] = 0.0
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),0.0)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
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
        pb.pbclass = "NLR2-AN-V-V"
        self.pb = pb; self.pbm = pbm

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

