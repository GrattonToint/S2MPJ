from s2mpjlib import *
class  DALLASS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DALLASS
#    *********
# 
#    The small Dallas water distribution problem
#    The problem is also named "W30" in some references.
#    This is a nonlinear network problem with conditioning of
#    the order of 10**4.
# 
#    Source:
#    R. Dembo,
#    private communication, 1986.
# 
#    SIF input: Ph. Toint, June 1990.
# 
#    classification = "C-CONR2-MN-46-31"
# 
#    Number of arcs
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'DALLASS'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 46
        v_['NODES'] = 31
        v_['1'] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X42']])
        valA = np.append(valA,float(-6.38400e+02))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X43']])
        valA = np.append(valA,float(-6.33000e+02))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X44']])
        valA = np.append(valA,float(-5.54500e+02))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X45']])
        valA = np.append(valA,float(-5.05000e+02))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X46']])
        valA = np.append(valA,float(-4.36900e+02))
        [ig,ig_,_] = s2mpj_ii('N1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N1')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X46']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X41']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X1']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('N2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N2')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X45']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('N3',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N3')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X44']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X3']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('N4',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N4')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X4']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('N5',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N5')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X16']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X7']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X6']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X5']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('N6',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N6')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X7']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X9']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X8']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('N7',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N7')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X9']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X10']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('N8',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N8')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X10']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X12']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X11']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('N9',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N9')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X12']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X3']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X13']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('N10',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N10')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X4']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X16']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X15']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X14']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('N11',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N11')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X15']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X13']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X5']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X17']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('N12',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N12')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X20']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X19']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X8']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X6']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X18']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('N13',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N13')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X42']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X18']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X19']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('N14',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N14')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X21']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X20']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('N15',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N15')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X43']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X21']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('N16',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N16')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X14']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X11']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X23']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X22']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('N17',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N17')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X23']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X25']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X24']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('N18',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N18')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X31']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X25']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X22']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X26']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('N19',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N19')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X26']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X17']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X28']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X27']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('N20',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N20')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X28']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('N21',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N21')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X31']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X30']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X29']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('N22',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N22')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X30']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X27']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('N23',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N23')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X24']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X32']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('N24',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N24')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X38']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X29']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X34']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X33']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('N25',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N25')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X32']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X35']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('N26',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N26')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X35']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X37']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X36']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('N27',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N27')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X37']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X34']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('N28',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N28')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X36']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X40']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X39']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X38']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('N29',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N29')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X39']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X33']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X1']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('N30',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N30')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X40']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X41']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('N31',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N31')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X46']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X45']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X44']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X43']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X42']])
        valA = np.append(valA,float(-1.0))
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
        self.gconst = arrset(self.gconst,ig_['N5'],float(2.80000))
        self.gconst = arrset(self.gconst,ig_['N7'],float(4.03000e-01))
        self.gconst = arrset(self.gconst,ig_['N8'],float(5.92000e-01))
        self.gconst = arrset(self.gconst,ig_['N9'],float(1.15600))
        self.gconst = arrset(self.gconst,ig_['N10'],float(2.00000e-01))
        self.gconst = arrset(self.gconst,ig_['N11'],float(4.95000e-01))
        self.gconst = arrset(self.gconst,ig_['N16'],float(3.13000e-01))
        self.gconst = arrset(self.gconst,ig_['N17'],float(8.44000e-01))
        self.gconst = arrset(self.gconst,ig_['N18'],float(3.31000e-01))
        self.gconst = arrset(self.gconst,ig_['N19'],float(5.30000e-02))
        self.gconst = arrset(self.gconst,ig_['N21'],float(2.72000e-01))
        self.gconst = arrset(self.gconst,ig_['N22'],float(8.83000e-01))
        self.gconst = arrset(self.gconst,ig_['N23'],float(5.71000e-01))
        self.gconst = arrset(self.gconst,ig_['N24'],float(7.55000e-01))
        self.gconst = arrset(self.gconst,ig_['N26'],float(5.27000e-01))
        self.gconst = arrset(self.gconst,ig_['N29'],float(1.00000e-03))
        self.gconst = arrset(self.gconst,ig_['N31'],float(-1.01960e+01))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-2.00000e+02)
        self.xupper = np.full((self.n,1),2.00000e+02)
        self.xlower[ix_['X1']] = 0.00000
        self.xupper[ix_['X1']] = 2.11673e+01
        self.xlower[ix_['X2']] = 0.00000
        self.xupper[ix_['X2']] = 4.37635e+01
        self.xlower[ix_['X3']] = 0.00000
        self.xupper[ix_['X3']] = 3.28255e+01
        self.xlower[ix_['X19']] = 0.00000
        self.xupper[ix_['X19']] = 2.20120e+01
        self.xlower[ix_['X21']] = 0.00000
        self.xupper[ix_['X21']] = 1.36703e+01
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(-2.00000e+02))
        if('X1' in ix_):
            self.x0[ix_['X1']] = float(2.11673e+01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X1']),float(2.11673e+01)))
        if('X2' in ix_):
            self.x0[ix_['X2']] = float(4.37635e+01)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X2'])[0],float(4.37635e+01)))
        if('X3' in ix_):
            self.x0[ix_['X3']] = float(3.28255e+01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X3']),float(3.28255e+01)))
        if('X4' in ix_):
            self.x0[ix_['X4']] = float(1.42109e-14)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X4'])[0],float(1.42109e-14)))
        if('X5' in ix_):
            self.x0[ix_['X5']] = float(1.68826e+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X5']),float(1.68826e+02)))
        if('X7' in ix_):
            self.x0[ix_['X7']] = float(2.81745e+01)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X7'])[0],float(2.81745e+01)))
        if('X8' in ix_):
            self.x0[ix_['X8']] = float(8.75603e+01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X8']),float(8.75603e+01)))
        if('X9' in ix_):
            self.x0[ix_['X9']] = float(-5.93858e+01)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X9'])[0],float(-5.93858e+01)))
        if('X10' in ix_):
            self.x0[ix_['X10']] = float(-5.97888e+01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X10']),float(-5.97888e+01)))
        if('X11' in ix_):
            self.x0[ix_['X11']] = float(1.83383e+02)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X11'])[0],float(1.83383e+02)))
        if('X13' in ix_):
            self.x0[ix_['X13']] = float(-1.68331e+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X13']),float(-1.68331e+02)))
        if('X15' in ix_):
            self.x0[ix_['X15']] = float(2.00000e+02)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X15'])[0],float(2.00000e+02)))
        if('X16' in ix_):
            self.x0[ix_['X16']] = float(2.00000e-01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X16']),float(2.00000e-01)))
        if('X17' in ix_):
            self.x0[ix_['X17']] = float(2.00000e+02)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X17'])[0],float(2.00000e+02)))
        if('X18' in ix_):
            self.x0[ix_['X18']] = float(-7.67574e+01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X18']),float(-7.67574e+01)))
        if('X19' in ix_):
            self.x0[ix_['X19']] = float(2.20120e+01)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X19'])[0],float(2.20120e+01)))
        if('X20' in ix_):
            self.x0[ix_['X20']] = float(1.36703e+01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X20']),float(1.36703e+01)))
        if('X21' in ix_):
            self.x0[ix_['X21']] = float(1.36703e+01)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X21'])[0],float(1.36703e+01)))
        if('X22' in ix_):
            self.x0[ix_['X22']] = float(-1.98461e+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X22']),float(-1.98461e+02)))
        if('X23' in ix_):
            self.x0[ix_['X23']] = float(1.81531e+02)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X23'])[0],float(1.81531e+02)))
        if('X24' in ix_):
            self.x0[ix_['X24']] = float(-1.93133e+01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X24']),float(-1.93133e+01)))
        if('X25' in ix_):
            self.x0[ix_['X25']] = float(2.00000e+02)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X25'])[0],float(2.00000e+02)))
        if('X26' in ix_):
            self.x0[ix_['X26']] = float(-1.98792e+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X26']),float(-1.98792e+02)))
        if('X27' in ix_):
            self.x0[ix_['X27']] = float(1.15500)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X27'])[0],float(1.15500)))
        if('X28' in ix_):
            self.x0[ix_['X28']] = float(0.00000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X28']),float(0.00000)))
        if('X29' in ix_):
            self.x0[ix_['X29']] = float(2.00000e+02)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X29'])[0],float(2.00000e+02)))
        if('X30' in ix_):
            self.x0[ix_['X30']] = float(2.72000e-01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X30']),float(2.72000e-01)))
        if('X32' in ix_):
            self.x0[ix_['X32']] = float(-1.98843e+01)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X32'])[0],float(-1.98843e+01)))
        if('X33' in ix_):
            self.x0[ix_['X33']] = float(1.78834e+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X33']),float(1.78834e+02)))
        if('X34' in ix_):
            self.x0[ix_['X34']] = float(-1.79589e+02)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X34'])[0],float(-1.79589e+02)))
        if('X35' in ix_):
            self.x0[ix_['X35']] = float(-1.98843e+01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X35']),float(-1.98843e+01)))
        if('X37' in ix_):
            self.x0[ix_['X37']] = float(1.79589e+02)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X37'])[0],float(1.79589e+02)))
        if('X40' in ix_):
            self.x0[ix_['X40']] = float(2.00000e+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X40']),float(2.00000e+02)))
        if('X41' in ix_):
            self.x0[ix_['X41']] = float(2.00000e+02)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X41'])[0],float(2.00000e+02)))
        if('X42' in ix_):
            self.x0[ix_['X42']] = float(9.87694e+01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X42']),float(9.87694e+01)))
        if('X43' in ix_):
            self.x0[ix_['X43']] = float(1.36703e+01)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X43'])[0],float(1.36703e+01)))
        if('X44' in ix_):
            self.x0[ix_['X44']] = float(3.28255e+01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X44']),float(3.28255e+01)))
        if('X45' in ix_):
            self.x0[ix_['X45']] = float(4.37635e+01)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X45'])[0],float(4.37635e+01)))
        if('X46' in ix_):
            self.x0[ix_['X46']] = float(-1.78833e+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X46']),float(-1.78833e+02)))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eT1', iet_)
        elftv = loaset(elftv,it,0,'ARC')
        elftp = []
        elftp = loaset(elftp,it,0,'C1')
        elftp = loaset(elftp,it,1,'C2')
        elftp = loaset(elftp,it,2,'C3')
        [it,iet_,_] = s2mpj_ii( 'eT4', iet_)
        elftv = loaset(elftv,it,0,'ARC')
        elftp = loaset(elftp,it,0,'C1')
        elftp = loaset(elftp,it,1,'C2')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        ename = 'E1'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eT4')
        ielftype = arrset(ielftype,ie,iet_["eT4"])
        vname = 'X1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(4.48060e+02))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(2.51200e+02))
        ename = 'E2'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eT4')
        ielftype = arrset(ielftype,ie,iet_["eT4"])
        vname = 'X2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.91526e+03))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(6.46300e+01))
        ename = 'E3'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eT4')
        ielftype = arrset(ielftype,ie,iet_["eT4"])
        vname = 'X3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.07752e+03))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(4.81400e+01))
        ename = 'E4'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eT1')
            ielftype = arrset(ielftype,ie,iet_['eT1'])
        vname = 'X4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(6.90000e+03))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(6.60000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.22000e+02))
        ename = 'E5'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eT1')
            ielftype = arrset(ielftype,ie,iet_['eT1'])
        vname = 'X5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.00000e+02))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.00000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.00000e+02))
        ename = 'E6'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eT1')
            ielftype = arrset(ielftype,ie,iet_['eT1'])
        vname = 'X6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(6.63000e+02))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(3.60000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.20000e+02))
        ename = 'E7'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eT1')
            ielftype = arrset(ielftype,ie,iet_['eT1'])
        vname = 'X7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(5.10000e+03))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(6.60000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.22000e+02))
        ename = 'E8'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eT1')
            ielftype = arrset(ielftype,ie,iet_['eT1'])
        vname = 'X8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(6.45000e+02))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(3.00000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.20000e+02))
        ename = 'E9'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eT1')
            ielftype = arrset(ielftype,ie,iet_['eT1'])
        vname = 'X9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(7.40000e+03))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(6.60000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.22000e+02))
        ename = 'E10'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eT1')
            ielftype = arrset(ielftype,ie,iet_['eT1'])
        vname = 'X10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(5.00000e+03))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(6.60000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(9.50000e+01))
        ename = 'E11'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eT1')
            ielftype = arrset(ielftype,ie,iet_['eT1'])
        vname = 'X11'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(8.00000e+02))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(5.40000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.07000e+02))
        ename = 'E12'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eT1')
            ielftype = arrset(ielftype,ie,iet_['eT1'])
        vname = 'X12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(5.20000e+03))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(4.80000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.10000e+02))
        ename = 'E13'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eT1')
            ielftype = arrset(ielftype,ie,iet_['eT1'])
        vname = 'X13'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(6.00000e+03))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(4.80000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.10000e+02))
        ename = 'E14'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eT1')
            ielftype = arrset(ielftype,ie,iet_['eT1'])
        vname = 'X14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(4.00000e+02))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(5.40000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.00000e+02))
        ename = 'E15'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eT1')
            ielftype = arrset(ielftype,ie,iet_['eT1'])
        vname = 'X15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(4.00000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(3.12200e+01))
        posep = np.where(elftp[ielftype[ie]]=='C3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.30000e+02))
        ename = 'E16'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eT1')
            ielftype = arrset(ielftype,ie,iet_['eT1'])
        vname = 'X16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(4.50000e+03))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(6.60000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.22000e+02))
        ename = 'E17'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eT1')
            ielftype = arrset(ielftype,ie,iet_['eT1'])
        vname = 'X17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(5.10000e+03))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(2.40000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.10000e+02))
        ename = 'E18'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eT1')
            ielftype = arrset(ielftype,ie,iet_['eT1'])
        vname = 'X18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(3.00000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(4.80000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.18000e+02))
        ename = 'E19'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eT4')
        ielftype = arrset(ielftype,ie,iet_["eT4"])
        vname = 'X19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(4.84530e+02))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.12970e+02))
        ename = 'E20'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eT1')
            ielftype = arrset(ielftype,ie,iet_['eT1'])
        vname = 'X20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(2.60000e+04))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(4.80000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.10000e+02))
        ename = 'E21'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eT4')
        ielftype = arrset(ielftype,ie,iet_["eT4"])
        vname = 'X21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.86880e+02))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.60610e+02))
        ename = 'E22'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eT1')
            ielftype = arrset(ielftype,ie,iet_['eT1'])
        vname = 'X22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(5.20000e+02))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(3.36100e+01))
        posep = np.where(elftp[ielftype[ie]]=='C3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.30000e+02))
        ename = 'E23'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eT1')
            ielftype = arrset(ielftype,ie,iet_['eT1'])
        vname = 'X23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(4.60000e+03))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(5.40000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(9.50000e+01))
        ename = 'E24'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eT1')
            ielftype = arrset(ielftype,ie,iet_['eT1'])
        vname = 'X24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(5.40000e+03))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(5.40000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(9.50000e+01))
        ename = 'E25'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eT1')
            ielftype = arrset(ielftype,ie,iet_['eT1'])
        vname = 'X25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(5.60000e+03))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.20000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.10000e+02))
        ename = 'E26'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eT1')
            ielftype = arrset(ielftype,ie,iet_['eT1'])
        vname = 'X26'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(3.30000e+03))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.20000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.10000e+02))
        ename = 'E27'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eT1')
            ielftype = arrset(ielftype,ie,iet_['eT1'])
        vname = 'X27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(2.20000e+03))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(2.40000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.24000e+02))
        ename = 'E28'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eT1')
            ielftype = arrset(ielftype,ie,iet_['eT1'])
        vname = 'X28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.00000e+03))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(2.40000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.10000e+02))
        ename = 'E29'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eT1')
            ielftype = arrset(ielftype,ie,iet_['eT1'])
        vname = 'X29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(5.90000e+03))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(2.40000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.13000e+02))
        ename = 'E30'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eT1')
            ielftype = arrset(ielftype,ie,iet_['eT1'])
        vname = 'X30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(2.80000e+03))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(2.40000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.13000e+02))
        ename = 'E31'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eT1')
            ielftype = arrset(ielftype,ie,iet_['eT1'])
        vname = 'X31'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(2.70000e+03))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.20000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.10000e+02))
        ename = 'E32'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eT1')
            ielftype = arrset(ielftype,ie,iet_['eT1'])
        vname = 'X32'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(3.10000e+03))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(5.40000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(9.50000e+01))
        ename = 'E33'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eT1')
            ielftype = arrset(ielftype,ie,iet_['eT1'])
        vname = 'X33'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.00000e+02))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.00000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.00000e+02))
        ename = 'E34'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eT1')
            ielftype = arrset(ielftype,ie,iet_['eT1'])
        vname = 'X34'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(4.30000e+03))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(2.40000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.13000e+02))
        ename = 'E35'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eT1')
            ielftype = arrset(ielftype,ie,iet_['eT1'])
        vname = 'X35'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(2.20000e+03))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(5.40000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(9.50000e+01))
        ename = 'E36'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eT1')
            ielftype = arrset(ielftype,ie,iet_['eT1'])
        vname = 'X36'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.80000e+03))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(5.40000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.10000e+02))
        ename = 'E37'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eT1')
            ielftype = arrset(ielftype,ie,iet_['eT1'])
        vname = 'X37'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.00000e+02))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(2.40000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.10000e+02))
        ename = 'E38'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eT1')
            ielftype = arrset(ielftype,ie,iet_['eT1'])
        vname = 'X38'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.31000e+03))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(3.00000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.00000e+02))
        ename = 'E39'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eT1')
            ielftype = arrset(ielftype,ie,iet_['eT1'])
        vname = 'X39'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(6.65000e+02))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(3.60000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.20000e+02))
        ename = 'E40'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eT1')
            ielftype = arrset(ielftype,ie,iet_['eT1'])
        vname = 'X40'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.10000e+03))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(3.60000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.20000e+02))
        ename = 'E41'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eT1')
            ielftype = arrset(ielftype,ie,iet_['eT1'])
        vname = 'X41'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-2.00000e+02),float(2.00000e+02),float(-2.00000e+02)))
        posev = np.where(elftv[ielftype[ie]]=='ARC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='C1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(3.23000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.00000e+01))
        posep = np.where(elftp[ielftype[ie]]=='C3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.00000e+02))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        ig = ig_['OBJ']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E1'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E2'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E3'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E4'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E5'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E6'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E7'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E8'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E9'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E10'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E11'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E12'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E13'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E14'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E15'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E16'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E17'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E18'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E19'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E20'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E21'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E22'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E23'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E24'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E25'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E26'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E27'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E28'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E29'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E30'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E31'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E32'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E33'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E34'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E35'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E36'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E37'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E38'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E39'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E40'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E41'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               -3.2393D+04
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CONR2-MN-46-31"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eT1(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        TMP = 850559.0e0/2.85*self.elpar[iel_][0]
        TMP = TMP/(self.elpar[iel_][2]**1.85)
        TMP = TMP/(self.elpar[iel_][1]**4.87)
        X = np.absolute(EV_[0])
        XEXP = X**0.85
        f_   = TMP*X**2*XEXP
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.85*TMP*EV_[0]*XEXP
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 5.2725*TMP*XEXP
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eT4(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        EPS2 = 1.0e-14
        SQC1 = np.sqrt(self.elpar[iel_][0])
        X = min(EV_[0],SQC1)
        TMP = self.elpar[iel_][1]*(self.elpar[iel_][0]-X*X)
        TMP = np.sqrt(max(TMP,EPS2))
        TMP2 = np.sqrt(self.elpar[iel_][1])*np.arcsin(X/SQC1)
        f_   = 0.5*(-X*TMP-self.elpar[iel_][0]*TMP2)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -TMP
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = self.elpar[iel_][1]*X/TMP
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

