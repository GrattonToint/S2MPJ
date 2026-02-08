from s2mpjlib import *
class  HIMMELBJ(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HIMMELBJ
#    *********
# 
#    An chemical equilibrium problem by A.P. Jones.
#    It has a nonlinear objective and linear constraints
# 
#    Source: problem 6 in
#    D.H. Himmelblau,
#    "Applied nonlinear programming",
#    McGraw-Hill, New-York, 1972.
# 
#    SIF input: Ph. Toint, March 1991.
# 
#    classification = "C-COLR2-MY-45-14"
# 
#    Number of variable sets
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HIMMELBJ'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['NSETS'] = 7
        v_['NS1'] = 4.0
        v_['NS2'] = 13.0
        v_['NS3'] = 18.0
        v_['NS4'] = 3.0
        v_['NS5'] = 3.0
        v_['NS6'] = 2.0
        v_['NS7'] = 2.0
        v_['NEQ'] = 16
        v_['1'] = 1
        v_['2'] = 2
        v_['3'] = 3
        v_['4'] = 4
        v_['13'] = 13
        v_['18'] = 18
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for K in range(int(v_['1']),int(v_['NSETS'])+1):
            v_['RNSK'] = v_['NS'+str(K)]
            v_['NSK'] = int(np.fix(v_['RNSK']))
            for J in range(int(v_['1']),int(v_['NSK'])+1):
                [iv,ix_,_] = s2mpj_ii('X'+str(J)+','+str(K),ix_)
                self.xnames=arrset(self.xnames,iv,'X'+str(J)+','+str(K))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2,1']])
        valA = np.append(valA,float(-7.69))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X3,1']])
        valA = np.append(valA,float(-11.52))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X4,1']])
        valA = np.append(valA,float(-36.60))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X1,2']])
        valA = np.append(valA,float(-10.94))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X8,2']])
        valA = np.append(valA,float(2.5966))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X9,2']])
        valA = np.append(valA,float(-39.39))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X10,2']])
        valA = np.append(valA,float(-21.35))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X11,2']])
        valA = np.append(valA,float(-32.84))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X12,2']])
        valA = np.append(valA,float(6.26))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X1,3']])
        valA = np.append(valA,float(10.45))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X3,3']])
        valA = np.append(valA,float(-0.5))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X7,3']])
        valA = np.append(valA,float(2.2435))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X9,3']])
        valA = np.append(valA,float(-39.39))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X10,3']])
        valA = np.append(valA,float(-21.49))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X11,3']])
        valA = np.append(valA,float(-32.84))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X12,3']])
        valA = np.append(valA,float(6.12))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X15,3']])
        valA = np.append(valA,float(-1.9028))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X16,3']])
        valA = np.append(valA,float(-2.8889))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X17,3']])
        valA = np.append(valA,float(-3.3622))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X18,3']])
        valA = np.append(valA,float(-7.4854))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X1,4']])
        valA = np.append(valA,float(-15.639))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X3,4']])
        valA = np.append(valA,float(21.81))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X1,5']])
        valA = np.append(valA,float(-16.79))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X3,5']])
        valA = np.append(valA,float(18.9779))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2,6']])
        valA = np.append(valA,float(11.959))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2,7']])
        valA = np.append(valA,float(12.899))
        [ig,ig_,_] = s2mpj_ii('C1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C1')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X1,1']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X1,2']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X1,3']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X15,3']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X16,3']])
        valA = np.append(valA,float(2.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X17,3']])
        valA = np.append(valA,float(3.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X18,3']])
        valA = np.append(valA,float(4.0))
        [ig,ig_,_] = s2mpj_ii('C2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C2')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2,1']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2,2']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X10,2']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X11,2']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X12,2']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2,3']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X10,3']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X11,3']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X12,3']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2,6']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2,7']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('C3',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C3')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X3,1']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X3,2']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X3,3']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('C4',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C4')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X4,1']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X4,2']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X9,2']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X11,2']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X12,2']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X4,3']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X9,3']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X11,3']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X12,3']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X1,4']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X3,4']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X1,5']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X3,5']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2,6']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2,7']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('C5',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C5')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X4,1']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X5,2']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X9,2']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X10,2']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X11,2']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X12,2']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X5,3']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X9,3']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X10,3']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X11,3']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X12,3']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('C6',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C6')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X6,2']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X6,3']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('C7',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C7')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X7,2']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X7,3']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('C8',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C8')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X8,2']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X8,3']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('C11',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C11')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X14,3']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X15,3']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X16,3']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X17,3']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X18,3']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('C12',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C12')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X4,2']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X5,2']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X6,2']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X7,2']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X8,2']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X10,2']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X12,2']])
        valA = np.append(valA,float(-2.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X13,2']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X14,3']])
        valA = np.append(valA,float(-4.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X15,3']])
        valA = np.append(valA,float(-3.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X16,3']])
        valA = np.append(valA,float(-2.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X17,3']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('C13',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C13')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X15,3']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X16,3']])
        valA = np.append(valA,float(-2.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X17,3']])
        valA = np.append(valA,float(-3.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X18,3']])
        valA = np.append(valA,float(-4.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X1,4']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2,4']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X3,4']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('C14',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C14')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X1,5']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2,5']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X3,5']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('C15',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C15')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X3,4']])
        valA = np.append(valA,float(-4.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X1,6']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2,6']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('C16',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C16')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X3,5']])
        valA = np.append(valA,float(-4.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X1,7']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2,7']])
        valA = np.append(valA,float(1.0))
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
        self.gconst = arrset(self.gconst,ig_['C1'],float(0.652981))
        self.gconst = arrset(self.gconst,ig_['C2'],float(0.281941))
        self.gconst = arrset(self.gconst,ig_['C3'],float(3.705233))
        self.gconst = arrset(self.gconst,ig_['C4'],float(47.00022))
        self.gconst = arrset(self.gconst,ig_['C5'],float(47.02972))
        self.gconst = arrset(self.gconst,ig_['C6'],float(0.08005))
        self.gconst = arrset(self.gconst,ig_['C7'],float(0.08813))
        self.gconst = arrset(self.gconst,ig_['C8'],float(0.04829))
        self.gconst = arrset(self.gconst,ig_['C11'],float(0.0022725))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),1.e-12)
        self.xupper = np.full((self.n,1),+float('inf'))
        self.xlower[ix_['X13,2']] = 0.0155
        self.xupper[ix_['X13,2']] = 0.0155
        self.xlower[ix_['X13,3']] = 0.0211275
        self.xupper[ix_['X13,3']] = 0.0211275
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(0.1))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eXLOGX', iet_)
        elftv = loaset(elftv,it,0,'X')
        [it,iet_,_] = s2mpj_ii( 'eXLOGX2', iet_)
        elftv = loaset(elftv,it,0,'Y1')
        elftv = loaset(elftv,it,1,'Y2')
        elftv = loaset(elftv,it,2,'X')
        [it,iet_,_] = s2mpj_ii( 'eXLOGX3', iet_)
        elftv = loaset(elftv,it,0,'Y1')
        elftv = loaset(elftv,it,1,'Y2')
        elftv = loaset(elftv,it,2,'Y3')
        elftv = loaset(elftv,it,3,'X')
        [it,iet_,_] = s2mpj_ii( 'eXLOGX4', iet_)
        elftv = loaset(elftv,it,0,'Y1')
        elftv = loaset(elftv,it,1,'Y2')
        elftv = loaset(elftv,it,2,'Y3')
        elftv = loaset(elftv,it,3,'Y4')
        elftv = loaset(elftv,it,4,'X')
        [it,iet_,_] = s2mpj_ii( 'eXLOGX13', iet_)
        elftv = loaset(elftv,it,0,'Y1')
        elftv = loaset(elftv,it,1,'Y2')
        elftv = loaset(elftv,it,2,'Y3')
        elftv = loaset(elftv,it,3,'Y4')
        elftv = loaset(elftv,it,4,'Y5')
        elftv = loaset(elftv,it,5,'Y6')
        elftv = loaset(elftv,it,6,'Y7')
        elftv = loaset(elftv,it,7,'Y8')
        elftv = loaset(elftv,it,8,'Y9')
        elftv = loaset(elftv,it,9,'Y10')
        elftv = loaset(elftv,it,10,'Y11')
        elftv = loaset(elftv,it,11,'Y12')
        elftv = loaset(elftv,it,12,'Y13')
        elftv = loaset(elftv,it,13,'X')
        [it,iet_,_] = s2mpj_ii( 'eXLOGX18', iet_)
        elftv = loaset(elftv,it,0,'Y1')
        elftv = loaset(elftv,it,1,'Y2')
        elftv = loaset(elftv,it,2,'Y3')
        elftv = loaset(elftv,it,3,'Y4')
        elftv = loaset(elftv,it,4,'Y5')
        elftv = loaset(elftv,it,5,'Y6')
        elftv = loaset(elftv,it,6,'Y7')
        elftv = loaset(elftv,it,7,'Y8')
        elftv = loaset(elftv,it,8,'Y9')
        elftv = loaset(elftv,it,9,'Y10')
        elftv = loaset(elftv,it,10,'Y11')
        elftv = loaset(elftv,it,11,'Y12')
        elftv = loaset(elftv,it,12,'Y13')
        elftv = loaset(elftv,it,13,'Y14')
        elftv = loaset(elftv,it,14,'Y15')
        elftv = loaset(elftv,it,15,'Y16')
        elftv = loaset(elftv,it,16,'Y17')
        elftv = loaset(elftv,it,17,'Y18')
        elftv = loaset(elftv,it,18,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for K in range(int(v_['1']),int(v_['NSETS'])+1):
            v_['RNSK'] = v_['NS'+str(K)]
            v_['NSK'] = int(np.fix(v_['RNSK']))
            for J in range(int(v_['1']),int(v_['NSK'])+1):
                ename = 'A'+str(J)+','+str(K)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eXLOGX')
                ielftype = arrset(ielftype,ie,iet_["eXLOGX"])
                vname = 'X'+str(J)+','+str(K)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
                posev = np.where(elftv[ielftype[ie]]=='X')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
        v_['K'] = 1
        for J in range(int(v_['1']),int(v_['4'])+1):
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eXLOGX4')
            ielftype = arrset(ielftype,ie,iet_["eXLOGX4"])
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(J)+','+str(int(v_['K']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X1,1'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X2,1'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X3,1'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X4,1'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y4')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        v_['K'] = 2
        for J in range(int(v_['1']),int(v_['13'])+1):
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eXLOGX13')
            ielftype = arrset(ielftype,ie,iet_["eXLOGX13"])
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(J)+','+str(int(v_['K']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X1,2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X2,2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X3,2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X4,2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y4')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X5,2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y5')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X6,2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y6')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X7,2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y7')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X8,2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y8')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X9,2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y9')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X10,2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y10')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X11,2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y11')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X12,2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y12')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X13,2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y13')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        v_['K'] = 3
        for J in range(int(v_['1']),int(v_['18'])+1):
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eXLOGX18')
            ielftype = arrset(ielftype,ie,iet_["eXLOGX18"])
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(J)+','+str(int(v_['K']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X1,3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X2,3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X3,3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X4,3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y4')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X5,3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y5')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X6,3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y6')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X7,3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y7')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X8,3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y8')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X9,3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y9')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X10,3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y10')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X11,3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y11')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X12,3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y12')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X13,3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y13')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X14,3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y14')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X15,3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y15')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X16,3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y16')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X17,3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y17')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X18,3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y18')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        v_['K'] = 4
        for J in range(int(v_['1']),int(v_['3'])+1):
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eXLOGX3')
            ielftype = arrset(ielftype,ie,iet_["eXLOGX3"])
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(J)+','+str(int(v_['K']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X1,4'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X2,4'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X3,4'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        v_['K'] = 5
        for J in range(int(v_['1']),int(v_['3'])+1):
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eXLOGX3')
            ielftype = arrset(ielftype,ie,iet_["eXLOGX3"])
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(J)+','+str(int(v_['K']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X1,5'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X2,5'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X3,5'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        v_['K'] = 6
        for J in range(int(v_['1']),int(v_['2'])+1):
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eXLOGX2')
            ielftype = arrset(ielftype,ie,iet_["eXLOGX2"])
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(J)+','+str(int(v_['K']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X1,6'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X2,6'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        v_['K'] = 7
        for J in range(int(v_['1']),int(v_['2'])+1):
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eXLOGX2')
            ielftype = arrset(ielftype,ie,iet_["eXLOGX2"])
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(J)+','+str(int(v_['K']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X1,7'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X2,7'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(1.e-12),None,float(0.1))
            posev = np.where(elftv[ielftype[ie]]=='Y2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for K in range(int(v_['1']),int(v_['NSETS'])+1):
            v_['RNSK'] = v_['NS'+str(K)]
            v_['NSK'] = int(np.fix(v_['RNSK']))
            for J in range(int(v_['1']),int(v_['NSK'])+1):
                ig = ig_['OBJ']
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['A'+str(J)+','+str(K)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,1.)
                posel = posel+1
                self.grelt = loaset(self.grelt,ig,posel,ie_['B'+str(J)+','+str(K)])
                self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               -1910.344724
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
        self.pbclass   = "C-COLR2-MY-45-14"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eXLOGX(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        LOGX = np.log(EV_[0,0])
        f_   = EV_[0,0]*LOGX
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = LOGX+1.0
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 1.0/EV_[0,0]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eXLOGX2(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((2,3))
        IV_ = np.zeros(2)
        U_[1,0] = U_[1,0]+1
        U_[1,1] = U_[1,1]+1
        U_[0,2] = U_[0,2]+1
        IV_[0] = to_scalar(U_[0:1,:].dot(EV_))
        IV_[1] = to_scalar(U_[1:2,:].dot(EV_))
        LOGX = np.log(IV_[1])
        f_   = IV_[0]*LOGX
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = LOGX
            g_[1] = IV_[0]/IV_[1]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = 1.0/IV_[1]
                H_[1,0] = H_[0,1]
                H_[1,1] = -IV_[0]/IV_[1]**2
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eXLOGX3(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((2,4))
        IV_ = np.zeros(2)
        U_[1,0] = U_[1,0]+1
        U_[1,1] = U_[1,1]+1
        U_[1,2] = U_[1,2]+1
        U_[0,3] = U_[0,3]+1
        IV_[0] = to_scalar(U_[0:1,:].dot(EV_))
        IV_[1] = to_scalar(U_[1:2,:].dot(EV_))
        LOGX = np.log(IV_[1])
        f_   = IV_[0]*LOGX
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = LOGX
            g_[1] = IV_[0]/IV_[1]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = 1.0/IV_[1]
                H_[1,0] = H_[0,1]
                H_[1,1] = -IV_[0]/IV_[1]**2
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eXLOGX4(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((2,5))
        IV_ = np.zeros(2)
        U_[1,0] = U_[1,0]+1
        U_[1,1] = U_[1,1]+1
        U_[1,2] = U_[1,2]+1
        U_[1,3] = U_[1,3]+1
        U_[0,4] = U_[0,4]+1
        IV_[0] = to_scalar(U_[0:1,:].dot(EV_))
        IV_[1] = to_scalar(U_[1:2,:].dot(EV_))
        LOGX = np.log(IV_[1])
        f_   = IV_[0]*LOGX
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = LOGX
            g_[1] = IV_[0]/IV_[1]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = 1.0/IV_[1]
                H_[1,0] = H_[0,1]
                H_[1,1] = -IV_[0]/IV_[1]**2
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eXLOGX13(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((2,14))
        IV_ = np.zeros(2)
        U_[1,0] = U_[1,0]+1
        U_[1,1] = U_[1,1]+1
        U_[1,2] = U_[1,2]+1
        U_[1,3] = U_[1,3]+1
        U_[1,4] = U_[1,4]+1
        U_[1,5] = U_[1,5]+1
        U_[1,6] = U_[1,6]+1
        U_[1,7] = U_[1,7]+1
        U_[1,8] = U_[1,8]+1
        U_[1,9] = U_[1,9]+1
        U_[1,10] = U_[1,10]+1
        U_[1,11] = U_[1,11]+1
        U_[1,12] = U_[1,12]+1
        U_[0,13] = U_[0,13]+1
        IV_[0] = to_scalar(U_[0:1,:].dot(EV_))
        IV_[1] = to_scalar(U_[1:2,:].dot(EV_))
        LOGX = np.log(IV_[1])
        f_   = IV_[0]*LOGX
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = LOGX
            g_[1] = IV_[0]/IV_[1]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = 1.0/IV_[1]
                H_[1,0] = H_[0,1]
                H_[1,1] = -IV_[0]/IV_[1]**2
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eXLOGX18(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((2,19))
        IV_ = np.zeros(2)
        U_[1,0] = U_[1,0]+1
        U_[1,1] = U_[1,1]+1
        U_[1,2] = U_[1,2]+1
        U_[1,3] = U_[1,3]+1
        U_[1,4] = U_[1,4]+1
        U_[1,5] = U_[1,5]+1
        U_[1,6] = U_[1,6]+1
        U_[1,7] = U_[1,7]+1
        U_[1,8] = U_[1,8]+1
        U_[1,9] = U_[1,9]+1
        U_[1,10] = U_[1,10]+1
        U_[1,11] = U_[1,11]+1
        U_[1,12] = U_[1,12]+1
        U_[1,13] = U_[1,13]+1
        U_[1,14] = U_[1,14]+1
        U_[1,15] = U_[1,15]+1
        U_[1,16] = U_[1,16]+1
        U_[1,17] = U_[1,17]+1
        U_[0,18] = U_[0,18]+1
        IV_[0] = to_scalar(U_[0:1,:].dot(EV_))
        IV_[1] = to_scalar(U_[1:2,:].dot(EV_))
        LOGX = np.log(IV_[1])
        f_   = IV_[0]*LOGX
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = LOGX
            g_[1] = IV_[0]/IV_[1]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = 1.0/IV_[1]
                H_[1,0] = H_[0,1]
                H_[1,1] = -IV_[0]/IV_[1]**2
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

