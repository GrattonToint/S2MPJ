from s2mpjlib import *
class  DEGENLPB(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DEGENLPB
#    *********
# 
#    A linear program with some degeneracy.
# 
#    Source:
#    T.C.T. Kotiah and D.I. Steinberg,
#    "Occurences of cycling and other phenomena arising in a class of
#    linear programming models",
#    Communications of the ACM, vol. 20, pp. 107-112, 1977.
# 
#    SIF input: Ph. Toint, Aug 1990.
# 
#    classification = "LLR2-AN-20-15"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'DEGENLPB'

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
        v_['N'] = 20
        v_['M'] = 15
        v_['1'] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X2']
        pbm.A[ig,iv] = float(-0.01)+pbm.A[ig,iv]
        iv = ix_['X3']
        pbm.A[ig,iv] = float(-33.333)+pbm.A[ig,iv]
        iv = ix_['X4']
        pbm.A[ig,iv] = float(-100.0)+pbm.A[ig,iv]
        iv = ix_['X5']
        pbm.A[ig,iv] = float(-0.01)+pbm.A[ig,iv]
        iv = ix_['X6']
        pbm.A[ig,iv] = float(-33.343)+pbm.A[ig,iv]
        iv = ix_['X7']
        pbm.A[ig,iv] = float(-100.01)+pbm.A[ig,iv]
        iv = ix_['X8']
        pbm.A[ig,iv] = float(-33.333)+pbm.A[ig,iv]
        iv = ix_['X9']
        pbm.A[ig,iv] = float(-133.33)+pbm.A[ig,iv]
        iv = ix_['X10']
        pbm.A[ig,iv] = float(-100.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C1')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X2']
        pbm.A[ig,iv] = float(2.0)+pbm.A[ig,iv]
        iv = ix_['X3']
        pbm.A[ig,iv] = float(2.0)+pbm.A[ig,iv]
        iv = ix_['X4']
        pbm.A[ig,iv] = float(2.0)+pbm.A[ig,iv]
        iv = ix_['X5']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X6']
        pbm.A[ig,iv] = float(2.0)+pbm.A[ig,iv]
        iv = ix_['X7']
        pbm.A[ig,iv] = float(2.0)+pbm.A[ig,iv]
        iv = ix_['X8']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X9']
        pbm.A[ig,iv] = float(2.0)+pbm.A[ig,iv]
        iv = ix_['X10']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C2')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X2']
        pbm.A[ig,iv] = float(300.0)+pbm.A[ig,iv]
        iv = ix_['X3']
        pbm.A[ig,iv] = float(0.09)+pbm.A[ig,iv]
        iv = ix_['X4']
        pbm.A[ig,iv] = float(0.03)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C3',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C3')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(0.326)+pbm.A[ig,iv]
        iv = ix_['X2']
        pbm.A[ig,iv] = float(-101.0)+pbm.A[ig,iv]
        iv = ix_['X5']
        pbm.A[ig,iv] = float(200.0)+pbm.A[ig,iv]
        iv = ix_['X6']
        pbm.A[ig,iv] = float(0.06)+pbm.A[ig,iv]
        iv = ix_['X7']
        pbm.A[ig,iv] = float(0.02)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C4',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C4')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(0.0066667)+pbm.A[ig,iv]
        iv = ix_['X3']
        pbm.A[ig,iv] = float(-1.03)+pbm.A[ig,iv]
        iv = ix_['X6']
        pbm.A[ig,iv] = float(200.0)+pbm.A[ig,iv]
        iv = ix_['X8']
        pbm.A[ig,iv] = float(0.06)+pbm.A[ig,iv]
        iv = ix_['X9']
        pbm.A[ig,iv] = float(0.02)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C5',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C5')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(6.6667e-4)+pbm.A[ig,iv]
        iv = ix_['X4']
        pbm.A[ig,iv] = float(-1.01)+pbm.A[ig,iv]
        iv = ix_['X7']
        pbm.A[ig,iv] = float(200.0)+pbm.A[ig,iv]
        iv = ix_['X9']
        pbm.A[ig,iv] = float(0.06)+pbm.A[ig,iv]
        iv = ix_['X10']
        pbm.A[ig,iv] = float(0.02)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C6',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C6')
        iv = ix_['X2']
        pbm.A[ig,iv] = float(0.978)+pbm.A[ig,iv]
        iv = ix_['X5']
        pbm.A[ig,iv] = float(-201.0)+pbm.A[ig,iv]
        iv = ix_['X11']
        pbm.A[ig,iv] = float(100.0)+pbm.A[ig,iv]
        iv = ix_['X12']
        pbm.A[ig,iv] = float(0.03)+pbm.A[ig,iv]
        iv = ix_['X13']
        pbm.A[ig,iv] = float(0.01)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C7',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C7')
        iv = ix_['X2']
        pbm.A[ig,iv] = float(0.01)+pbm.A[ig,iv]
        iv = ix_['X3']
        pbm.A[ig,iv] = float(0.489)+pbm.A[ig,iv]
        iv = ix_['X6']
        pbm.A[ig,iv] = float(-101.03)+pbm.A[ig,iv]
        iv = ix_['X12']
        pbm.A[ig,iv] = float(100.0)+pbm.A[ig,iv]
        iv = ix_['X14']
        pbm.A[ig,iv] = float(0.03)+pbm.A[ig,iv]
        iv = ix_['X15']
        pbm.A[ig,iv] = float(0.01)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C8',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C8')
        iv = ix_['X2']
        pbm.A[ig,iv] = float(0.001)+pbm.A[ig,iv]
        iv = ix_['X4']
        pbm.A[ig,iv] = float(0.489)+pbm.A[ig,iv]
        iv = ix_['X7']
        pbm.A[ig,iv] = float(-101.03)+pbm.A[ig,iv]
        iv = ix_['X13']
        pbm.A[ig,iv] = float(100.0)+pbm.A[ig,iv]
        iv = ix_['X15']
        pbm.A[ig,iv] = float(0.03)+pbm.A[ig,iv]
        iv = ix_['X16']
        pbm.A[ig,iv] = float(0.01)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C9',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C9')
        iv = ix_['X3']
        pbm.A[ig,iv] = float(0.001)+pbm.A[ig,iv]
        iv = ix_['X4']
        pbm.A[ig,iv] = float(0.01)+pbm.A[ig,iv]
        iv = ix_['X9']
        pbm.A[ig,iv] = float(-1.04)+pbm.A[ig,iv]
        iv = ix_['X15']
        pbm.A[ig,iv] = float(100.0)+pbm.A[ig,iv]
        iv = ix_['X18']
        pbm.A[ig,iv] = float(0.03)+pbm.A[ig,iv]
        iv = ix_['X19']
        pbm.A[ig,iv] = float(0.01)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C10',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C10')
        iv = ix_['X3']
        pbm.A[ig,iv] = float(0.02)+pbm.A[ig,iv]
        iv = ix_['X8']
        pbm.A[ig,iv] = float(-1.06)+pbm.A[ig,iv]
        iv = ix_['X14']
        pbm.A[ig,iv] = float(100.0)+pbm.A[ig,iv]
        iv = ix_['X17']
        pbm.A[ig,iv] = float(0.03)+pbm.A[ig,iv]
        iv = ix_['X19']
        pbm.A[ig,iv] = float(0.01)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C11',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C11')
        iv = ix_['X4']
        pbm.A[ig,iv] = float(0.002)+pbm.A[ig,iv]
        iv = ix_['X10']
        pbm.A[ig,iv] = float(-1.02)+pbm.A[ig,iv]
        iv = ix_['X16']
        pbm.A[ig,iv] = float(100.0)+pbm.A[ig,iv]
        iv = ix_['X19']
        pbm.A[ig,iv] = float(0.03)+pbm.A[ig,iv]
        iv = ix_['X20']
        pbm.A[ig,iv] = float(0.01)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C12',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C12')
        iv = ix_['X11']
        pbm.A[ig,iv] = float(-2.5742e-6)+pbm.A[ig,iv]
        iv = ix_['X13']
        pbm.A[ig,iv] = float(0.00252)+pbm.A[ig,iv]
        iv = ix_['X16']
        pbm.A[ig,iv] = float(-0.61975)+pbm.A[ig,iv]
        iv = ix_['X20']
        pbm.A[ig,iv] = float(1.03)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C13',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C13')
        iv = ix_['X11']
        pbm.A[ig,iv] = float(-0.00257)+pbm.A[ig,iv]
        iv = ix_['X12']
        pbm.A[ig,iv] = float(0.25221)+pbm.A[ig,iv]
        iv = ix_['X14']
        pbm.A[ig,iv] = float(-6.2)+pbm.A[ig,iv]
        iv = ix_['X17']
        pbm.A[ig,iv] = float(1.09)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C14',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C14')
        iv = ix_['X11']
        pbm.A[ig,iv] = float(0.00629)+pbm.A[ig,iv]
        iv = ix_['X12']
        pbm.A[ig,iv] = float(-0.20555)+pbm.A[ig,iv]
        iv = ix_['X13']
        pbm.A[ig,iv] = float(-4.1106)+pbm.A[ig,iv]
        iv = ix_['X15']
        pbm.A[ig,iv] = float(101.04)+pbm.A[ig,iv]
        iv = ix_['X16']
        pbm.A[ig,iv] = float(505.1)+pbm.A[ig,iv]
        iv = ix_['X19']
        pbm.A[ig,iv] = float(-256.72)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C15',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C15')
        iv = ix_['X12']
        pbm.A[ig,iv] = float(0.00841)+pbm.A[ig,iv]
        iv = ix_['X13']
        pbm.A[ig,iv] = float(-0.08406)+pbm.A[ig,iv]
        iv = ix_['X14']
        pbm.A[ig,iv] = float(-0.20667)+pbm.A[ig,iv]
        iv = ix_['X16']
        pbm.A[ig,iv] = float(20.658)+pbm.A[ig,iv]
        iv = ix_['X18']
        pbm.A[ig,iv] = float(1.07)+pbm.A[ig,iv]
        iv = ix_['X19']
        pbm.A[ig,iv] = float(-10.5)+pbm.A[ig,iv]
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
        pbm.gconst = arrset(pbm.gconst,ig_['C1'],float(0.70785))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xupper = np.full((pb.n,1),1.0)
        pb.xlower = np.zeros((pb.n,1))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN               3.06435
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
        pb.pbclass = "LLR2-AN-20-15"
        self.pb = pb; self.pbm = pbm

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

