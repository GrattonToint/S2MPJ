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
#    classification = "C-CLLR2-AN-20-15"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 17 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'DEGENLPB'

    def __init__(self, *args): 
        import numpy as np
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 20
        v_['M'] = 15
        v_['1'] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
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
        iv = ix_['X2']
        self.A[ig,iv] = float(-0.01)+self.A[ig,iv]
        iv = ix_['X3']
        self.A[ig,iv] = float(-33.333)+self.A[ig,iv]
        iv = ix_['X4']
        self.A[ig,iv] = float(-100.0)+self.A[ig,iv]
        iv = ix_['X5']
        self.A[ig,iv] = float(-0.01)+self.A[ig,iv]
        iv = ix_['X6']
        self.A[ig,iv] = float(-33.343)+self.A[ig,iv]
        iv = ix_['X7']
        self.A[ig,iv] = float(-100.01)+self.A[ig,iv]
        iv = ix_['X8']
        self.A[ig,iv] = float(-33.333)+self.A[ig,iv]
        iv = ix_['X9']
        self.A[ig,iv] = float(-133.33)+self.A[ig,iv]
        iv = ix_['X10']
        self.A[ig,iv] = float(-100.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C1')
        iv = ix_['X1']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X2']
        self.A[ig,iv] = float(2.0)+self.A[ig,iv]
        iv = ix_['X3']
        self.A[ig,iv] = float(2.0)+self.A[ig,iv]
        iv = ix_['X4']
        self.A[ig,iv] = float(2.0)+self.A[ig,iv]
        iv = ix_['X5']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X6']
        self.A[ig,iv] = float(2.0)+self.A[ig,iv]
        iv = ix_['X7']
        self.A[ig,iv] = float(2.0)+self.A[ig,iv]
        iv = ix_['X8']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X9']
        self.A[ig,iv] = float(2.0)+self.A[ig,iv]
        iv = ix_['X10']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C2')
        iv = ix_['X1']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X2']
        self.A[ig,iv] = float(300.0)+self.A[ig,iv]
        iv = ix_['X3']
        self.A[ig,iv] = float(0.09)+self.A[ig,iv]
        iv = ix_['X4']
        self.A[ig,iv] = float(0.03)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C3',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C3')
        iv = ix_['X1']
        self.A[ig,iv] = float(0.326)+self.A[ig,iv]
        iv = ix_['X2']
        self.A[ig,iv] = float(-101.0)+self.A[ig,iv]
        iv = ix_['X5']
        self.A[ig,iv] = float(200.0)+self.A[ig,iv]
        iv = ix_['X6']
        self.A[ig,iv] = float(0.06)+self.A[ig,iv]
        iv = ix_['X7']
        self.A[ig,iv] = float(0.02)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C4',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C4')
        iv = ix_['X1']
        self.A[ig,iv] = float(0.0066667)+self.A[ig,iv]
        iv = ix_['X3']
        self.A[ig,iv] = float(-1.03)+self.A[ig,iv]
        iv = ix_['X6']
        self.A[ig,iv] = float(200.0)+self.A[ig,iv]
        iv = ix_['X8']
        self.A[ig,iv] = float(0.06)+self.A[ig,iv]
        iv = ix_['X9']
        self.A[ig,iv] = float(0.02)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C5',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C5')
        iv = ix_['X1']
        self.A[ig,iv] = float(6.6667e-4)+self.A[ig,iv]
        iv = ix_['X4']
        self.A[ig,iv] = float(-1.01)+self.A[ig,iv]
        iv = ix_['X7']
        self.A[ig,iv] = float(200.0)+self.A[ig,iv]
        iv = ix_['X9']
        self.A[ig,iv] = float(0.06)+self.A[ig,iv]
        iv = ix_['X10']
        self.A[ig,iv] = float(0.02)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C6',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C6')
        iv = ix_['X2']
        self.A[ig,iv] = float(0.978)+self.A[ig,iv]
        iv = ix_['X5']
        self.A[ig,iv] = float(-201.0)+self.A[ig,iv]
        iv = ix_['X11']
        self.A[ig,iv] = float(100.0)+self.A[ig,iv]
        iv = ix_['X12']
        self.A[ig,iv] = float(0.03)+self.A[ig,iv]
        iv = ix_['X13']
        self.A[ig,iv] = float(0.01)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C7',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C7')
        iv = ix_['X2']
        self.A[ig,iv] = float(0.01)+self.A[ig,iv]
        iv = ix_['X3']
        self.A[ig,iv] = float(0.489)+self.A[ig,iv]
        iv = ix_['X6']
        self.A[ig,iv] = float(-101.03)+self.A[ig,iv]
        iv = ix_['X12']
        self.A[ig,iv] = float(100.0)+self.A[ig,iv]
        iv = ix_['X14']
        self.A[ig,iv] = float(0.03)+self.A[ig,iv]
        iv = ix_['X15']
        self.A[ig,iv] = float(0.01)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C8',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C8')
        iv = ix_['X2']
        self.A[ig,iv] = float(0.001)+self.A[ig,iv]
        iv = ix_['X4']
        self.A[ig,iv] = float(0.489)+self.A[ig,iv]
        iv = ix_['X7']
        self.A[ig,iv] = float(-101.03)+self.A[ig,iv]
        iv = ix_['X13']
        self.A[ig,iv] = float(100.0)+self.A[ig,iv]
        iv = ix_['X15']
        self.A[ig,iv] = float(0.03)+self.A[ig,iv]
        iv = ix_['X16']
        self.A[ig,iv] = float(0.01)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C9',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C9')
        iv = ix_['X3']
        self.A[ig,iv] = float(0.001)+self.A[ig,iv]
        iv = ix_['X4']
        self.A[ig,iv] = float(0.01)+self.A[ig,iv]
        iv = ix_['X9']
        self.A[ig,iv] = float(-1.04)+self.A[ig,iv]
        iv = ix_['X15']
        self.A[ig,iv] = float(100.0)+self.A[ig,iv]
        iv = ix_['X18']
        self.A[ig,iv] = float(0.03)+self.A[ig,iv]
        iv = ix_['X19']
        self.A[ig,iv] = float(0.01)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C10',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C10')
        iv = ix_['X3']
        self.A[ig,iv] = float(0.02)+self.A[ig,iv]
        iv = ix_['X8']
        self.A[ig,iv] = float(-1.06)+self.A[ig,iv]
        iv = ix_['X14']
        self.A[ig,iv] = float(100.0)+self.A[ig,iv]
        iv = ix_['X17']
        self.A[ig,iv] = float(0.03)+self.A[ig,iv]
        iv = ix_['X19']
        self.A[ig,iv] = float(0.01)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C11',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C11')
        iv = ix_['X4']
        self.A[ig,iv] = float(0.002)+self.A[ig,iv]
        iv = ix_['X10']
        self.A[ig,iv] = float(-1.02)+self.A[ig,iv]
        iv = ix_['X16']
        self.A[ig,iv] = float(100.0)+self.A[ig,iv]
        iv = ix_['X19']
        self.A[ig,iv] = float(0.03)+self.A[ig,iv]
        iv = ix_['X20']
        self.A[ig,iv] = float(0.01)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C12',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C12')
        iv = ix_['X11']
        self.A[ig,iv] = float(-2.5742e-6)+self.A[ig,iv]
        iv = ix_['X13']
        self.A[ig,iv] = float(0.00252)+self.A[ig,iv]
        iv = ix_['X16']
        self.A[ig,iv] = float(-0.61975)+self.A[ig,iv]
        iv = ix_['X20']
        self.A[ig,iv] = float(1.03)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C13',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C13')
        iv = ix_['X11']
        self.A[ig,iv] = float(-0.00257)+self.A[ig,iv]
        iv = ix_['X12']
        self.A[ig,iv] = float(0.25221)+self.A[ig,iv]
        iv = ix_['X14']
        self.A[ig,iv] = float(-6.2)+self.A[ig,iv]
        iv = ix_['X17']
        self.A[ig,iv] = float(1.09)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C14',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C14')
        iv = ix_['X11']
        self.A[ig,iv] = float(0.00629)+self.A[ig,iv]
        iv = ix_['X12']
        self.A[ig,iv] = float(-0.20555)+self.A[ig,iv]
        iv = ix_['X13']
        self.A[ig,iv] = float(-4.1106)+self.A[ig,iv]
        iv = ix_['X15']
        self.A[ig,iv] = float(101.04)+self.A[ig,iv]
        iv = ix_['X16']
        self.A[ig,iv] = float(505.1)+self.A[ig,iv]
        iv = ix_['X19']
        self.A[ig,iv] = float(-256.72)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C15',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C15')
        iv = ix_['X12']
        self.A[ig,iv] = float(0.00841)+self.A[ig,iv]
        iv = ix_['X13']
        self.A[ig,iv] = float(-0.08406)+self.A[ig,iv]
        iv = ix_['X14']
        self.A[ig,iv] = float(-0.20667)+self.A[ig,iv]
        iv = ix_['X16']
        self.A[ig,iv] = float(20.658)+self.A[ig,iv]
        iv = ix_['X18']
        self.A[ig,iv] = float(1.07)+self.A[ig,iv]
        iv = ix_['X19']
        self.A[ig,iv] = float(-10.5)+self.A[ig,iv]
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
        self.cnames= cnames[self.congrps]
        self.nob = ngrp-self.m
        self.objgrps = np.where(gtype=='<>')[0]
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        self.gconst = arrset(self.gconst,ig_['C1'],float(0.70785))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xupper = np.full((self.n,1),1.0)
        self.xlower = np.zeros((self.n,1))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 0.0
#    Solution
# LO SOLTN               3.06435
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        self.A.resize(ngrp,self.n)
        self.A     = self.A.tocsr()
        sA1,sA2    = self.A.shape
        self.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons   = np.arange(len(self.congrps))
        self.pbclass = "C-CLLR2-AN-20-15"
        self.objderlvl = 2
        self.conderlvl = [2]


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

