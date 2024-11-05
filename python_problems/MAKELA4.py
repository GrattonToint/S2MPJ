from s2mpjlib import *
class  MAKELA4(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : MAKELA4
#    *********
# 
#    A nonlinear minmax problem in twenty variables.
# 
#    Source: 
#    M.M. Makela,
#    "Nonsmooth optimization",
#    Ph.D. thesis, Jyvaskyla University, 1990
# 
#    SIF input: Ph. Toint, Nov 1993.
# 
#    classification = "C-CLLR2-AN-21-40"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 17 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'MAKELA4'

    def __init__(self, *args): 
        import numpy as np
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['1'] = 1
        v_['20'] = 20
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['20'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I))
        [iv,ix_,_] = s2mpj_ii('U',ix_)
        self.xnames=arrset(self.xnames,iv,'U')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.A       = lil_matrix((1000000,1000000))
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames      = np.array([])
        self.cnames = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['U']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        for I in range(int(v_['1']),int(v_['20'])+1):
            [ig,ig_,_] = s2mpj_ii('F'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'F'+str(I))
            iv = ix_['U']
            self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
            iv = ix_['X'+str(I)]
            self.A[ig,iv] = float(1.0)+self.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('MF'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'MF'+str(I))
            iv = ix_['U']
            self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
            iv = ix_['X'+str(I)]
            self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
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
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        self.xlower = np.zeros((self.n,1))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        if('X1' in ix_):
            self.x0[ix_['X1']] = float(1.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X1']),float(1.0)))
        if('X2' in ix_):
            self.x0[ix_['X2']] = float(2.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X2']),float(2.0)))
        if('X3' in ix_):
            self.x0[ix_['X3']] = float(3.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X3']),float(3.0)))
        if('X4' in ix_):
            self.x0[ix_['X4']] = float(4.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X4']),float(4.0)))
        if('X5' in ix_):
            self.x0[ix_['X5']] = float(5.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X5']),float(5.0)))
        if('X6' in ix_):
            self.x0[ix_['X6']] = float(6.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X6']),float(6.0)))
        if('X7' in ix_):
            self.x0[ix_['X7']] = float(7.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X7']),float(7.0)))
        if('X8' in ix_):
            self.x0[ix_['X8']] = float(8.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X8']),float(8.0)))
        if('X9' in ix_):
            self.x0[ix_['X9']] = float(9.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X9']),float(9.0)))
        if('X10' in ix_):
            self.x0[ix_['X10']] = float(10.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X10']),float(10.0)))
        if('X11' in ix_):
            self.x0[ix_['X11']] = float(-11.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X11']),float(-11.0)))
        if('X12' in ix_):
            self.x0[ix_['X12']] = float(-12.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X12']),float(-12.0)))
        if('X13' in ix_):
            self.x0[ix_['X13']] = float(-13.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X13']),float(-13.0)))
        if('X14' in ix_):
            self.x0[ix_['X14']] = float(-14.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X14']),float(-14.0)))
        if('X15' in ix_):
            self.x0[ix_['X15']] = float(-15.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X15']),float(-15.0)))
        if('X16' in ix_):
            self.x0[ix_['X16']] = float(-16.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X16']),float(-16.0)))
        if('X17' in ix_):
            self.x0[ix_['X17']] = float(-17.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X17']),float(-17.0)))
        if('X18' in ix_):
            self.x0[ix_['X18']] = float(-18.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X18']),float(-18.0)))
        if('X19' in ix_):
            self.x0[ix_['X19']] = float(-19.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X19']),float(-19.0)))
        if('X20' in ix_):
            self.x0[ix_['X20']] = float(-20.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X20']),float(-20.0)))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.cupper[np.arange(self.nle)] = np.zeros((self.nle,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        self.A.resize(ngrp,self.n)
        self.A     = self.A.tocsr()
        sA1,sA2    = self.A.shape
        self.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons   = np.arange(len(self.congrps))
        self.pbclass = "C-CLLR2-AN-21-40"
        self.objderlvl = 2
        self.conderlvl = [2]


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

