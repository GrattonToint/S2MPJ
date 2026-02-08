from s2mpjlib import *
class  MOSARQP2(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : MOSARQP2
#    *********
# 
#    A convex quadratic problem, with variable dimensions.
#    In this problem, a third of the linear constraints are active at the
#    solution. 
# 
#    Source:
#    J.L. Morales-Perez and R.W.H. Sargent,
#    "On the implementation and performance of an interior point method for
#    large sparse convex quadratic programming",
#    Centre for Process Systems Engineering, Imperial College, London,
#    November 1991.
# 
#    SIF input: Ph. Toint, August 1993.
#               minor correction by Ph. Shott, Jan 1995.
# 
#    classification = "C-CQLR2-AN-V-V"
# 
#    Problem variants: these are distinguished by the triplet ( N, M, COND ),
#    where: - N (nb of variables) must be a multiple of 3
#             and have an integer square root
#           - M (nb of constraints) must be at least sqrt(N) 
#             and at most N - sqrt(N)
#           - COND (problem conditioning) is a positive integer
#    Except for the first, the instances suggested are those used by Morales
#    and Sargent.
# 
#           Alternative values for the SIF file parameters:
# IE N                   36             $-PARAMETER
# IE M                   10             $-PARAMETER
# RE COND                2.0            $-PARAMETER
# 
# IE N                   100            $-PARAMETER     original value
# IE M                   10             $-PARAMETER     original value
# RE COND                3.0            $-PARAMETER     original value
# 
# IE N                   900            $-PARAMETER
# IE M                   30             $-PARAMETER
# RE COND                1.0            $-PARAMETER
# 
# IE N                   900            $-PARAMETER
# IE M                   30             $-PARAMETER
# RE COND                2.0            $-PARAMETER
# 
# IE N                   900            $-PARAMETER
# IE M                   30             $-PARAMETER
# RE COND                3.0            $-PARAMETER
# 
# IE N                   900            $-PARAMETER
# IE M                   60             $-PARAMETER
# RE COND                1.0            $-PARAMETER
# 
# IE N                   900            $-PARAMETER
# IE M                   60             $-PARAMETER
# RE COND                2.0            $-PARAMETER
# 
# IE N                   900            $-PARAMETER
# IE M                   60             $-PARAMETER
# RE COND                3.0            $-PARAMETER
# 
# IE N                   900            $-PARAMETER
# IE M                   90             $-PARAMETER
# RE COND                1.0            $-PARAMETER
# 
# IE N                   900            $-PARAMETER
# IE M                   90             $-PARAMETER
# RE COND                2.0            $-PARAMETER
# 
# IE N                   900            $-PARAMETER
# IE M                   90             $-PARAMETER
# RE COND                3.0            $-PARAMETER
# 
# IE N                   900            $-PARAMETER
# IE M                   120            $-PARAMETER
# RE COND                1.0            $-PARAMETER
# 
# IE N                   900            $-PARAMETER
# IE M                   120            $-PARAMETER
# RE COND                2.0            $-PARAMETER
# 
# IE N                   900            $-PARAMETER
# IE M                   120            $-PARAMETER
# RE COND                3.0            $-PARAMETER
# 
# IE N                   900            $-PARAMETER
# IE M                   300            $-PARAMETER
# RE COND                1.0            $-PARAMETER
# 
# IE N                   900            $-PARAMETER
# IE M                   300            $-PARAMETER
# RE COND                2.0            $-PARAMETER
# 
# IE N                   900            $-PARAMETER
# IE M                   300            $-PARAMETER
# RE COND                3.0            $-PARAMETER
# 
# IE N                   900            $-PARAMETER
# IE M                   600            $-PARAMETER
# RE COND                1.0            $-PARAMETER
# 
# IE N                   900            $-PARAMETER
# IE M                   600            $-PARAMETER
# RE COND                2.0            $-PARAMETER
# 
# IE N                   900            $-PARAMETER
# IE M                   600            $-PARAMETER
# RE COND                3.0            $-PARAMETER
# 
# IE N                   2500           $-PARAMETER
# IE M                   700            $-PARAMETER
# RE COND                1.0            $-PARAMETER
# 
# IE N                   2500           $-PARAMETER
# IE M                   700            $-PARAMETER
# RE COND                2.0            $-PARAMETER
# 
# IE N                   2500           $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'MOSARQP2'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(36);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
# IE M                   700            $-PARAMETER
        if nargin<2:
            v_['M'] = int(10);  #  SIF file default value
        else:
            v_['M'] = int(args[1])
# RE COND                3.0            $-PARAMETER
        if nargin<3:
            v_['COND'] = float(2.0);  #  SIF file default value
        else:
            v_['COND'] = float(args[2])
        v_['1'] = 1
        v_['2'] = 2
        v_['3'] = 3
        v_['N-1'] = -1+v_['N']
        v_['N-2'] = -2+v_['N']
        v_['RN-1'] = float(v_['N-1'])
        v_['RN'] = float(v_['N'])
        v_['M-1'] = -1+v_['M']
        v_['RNP'] = 0.1+v_['RN']
        v_['RRTN'] = np.sqrt(v_['RNP'])
        v_['RTN'] = int(np.fix(v_['RRTN']))
        v_['RTN-1'] = -1+v_['RTN']
        v_['RTN-2'] = -2+v_['RTN']
        v_['RTN+1'] = 1+v_['RTN']
        v_['2RTN-1'] = v_['RTN']+v_['RTN-1']
        v_['2RTN'] = v_['RTN']+v_['RTN']
        v_['M-RTN+1'] = v_['M']-v_['RTN-1']
        for I in range(int(v_['1']),int(v_['N-2'])+1,int(v_['3'])):
            v_['I+1'] = 1+I
            v_['I+2'] = 2+I
            v_['XC'+str(I)] = -1.0
            v_['XC'+str(int(v_['I+1']))] = -1.0
            v_['XC'+str(int(v_['I+2']))] = 1.0
        v_['XC'+str(int(v_['N-1']))] = -1.0
        v_['XC'+str(int(v_['N']))] = 1.0
        v_['NNZ'] = 10
        v_['Y1'] = -0.3569732
        v_['Y2'] = 0.9871576
        v_['Y3'] = 0.5619363
        v_['Y4'] = -0.1984624
        v_['Y5'] = 0.4653328
        v_['Y6'] = 0.7364367
        v_['Y7'] = -0.4560378
        v_['Y8'] = -0.6457813
        v_['Y9'] = -0.0601357
        v_['Y10'] = 0.1035624
        v_['NZ1'] = 0.68971452
        v_['NZ2'] = 0.13452678
        v_['NZ3'] = 0.51234678
        v_['NZ4'] = 0.76591423
        v_['NZ5'] = 0.20857854
        v_['NZ6'] = 0.85672348
        v_['NZ7'] = 0.04356789
        v_['NZ8'] = 0.44692743
        v_['NZ9'] = 0.30136413
        v_['NZ10'] = 0.91367489
        v_['YN2'] = 0.0
        for I in range(int(v_['1']),int(v_['NNZ'])+1):
            v_['RKI'] = v_['NZ'+str(I)]*v_['RN']
            v_['K'+str(I)] = 1.1+v_['RKI']
            v_['TMP'] = v_['Y'+str(I)]*v_['Y'+str(I)]
            v_['YN2'] = v_['YN2']+v_['TMP']
        v_['-2/YN2'] = -2.0/v_['YN2']
        v_['4/YN4'] = v_['-2/YN2']*v_['-2/YN2']
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['I-1'] = -1+I
            v_['RI-1'] = float(v_['I-1'])
            v_['TMP'] = v_['RI-1']/v_['RN-1']
            v_['TMP'] = v_['TMP']*v_['COND']
            v_['D'+str(I)] = np.exp(v_['TMP'])
        v_['YDY'] = 0.0
        v_['YXC'] = 0.0
        v_['YDXC'] = 0.0
        for I in range(int(v_['1']),int(v_['NNZ'])+1):
            v_['RKI'] = v_['K'+str(I)]
            v_['KI'] = int(np.fix(v_['RKI']))
            v_['DY'+str(I)] = v_['Y'+str(I)]*v_['D'+str(int(v_['KI']))]
            v_['TMP'] = v_['DY'+str(I)]*v_['Y'+str(I)]
            v_['YDY'] = v_['YDY']+v_['TMP']
            v_['TMP'] = v_['Y'+str(I)]*v_['XC'+str(int(v_['KI']))]
            v_['YXC'] = v_['YXC']+v_['TMP']
            v_['TMP'] = v_['DY'+str(I)]*v_['XC'+str(int(v_['KI']))]
            v_['YDXC'] = v_['YDXC']+v_['TMP']
        v_['AA'] = v_['-2/YN2']*v_['YXC']
        v_['DD'] = v_['4/YN4']*v_['YDY']
        v_['BB'] = v_['DD']*v_['YXC']
        v_['CC'] = v_['-2/YN2']*v_['YDXC']
        v_['BB+CC'] = v_['BB']+v_['CC']
        v_['DD/2'] = 0.5*v_['DD']
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['C'+str(I)] = v_['D'+str(I)]*v_['XC'+str(I)]
        for I in range(int(v_['1']),int(v_['NNZ'])+1):
            v_['RKI'] = v_['K'+str(I)]
            v_['KI'] = int(np.fix(v_['RKI']))
            v_['TMP'] = v_['DY'+str(I)]*v_['AA']
            v_['C'+str(int(v_['KI']))] = v_['C'+str(int(v_['KI']))]+v_['TMP']
            v_['TMP'] = v_['Y'+str(I)]*v_['BB+CC']
            v_['C'+str(int(v_['KI']))] = v_['C'+str(int(v_['KI']))]+v_['TMP']
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
        for I in range(int(v_['1']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(I)]])
            valA = np.append(valA,float(v_['C'+str(I)]))
        [ig,ig_,_] = s2mpj_ii('CS'+str(int(v_['1'])),ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'CS'+str(int(v_['1'])))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X'+str(int(v_['1']))]])
        valA = np.append(valA,float(4.0))
        [ig,ig_,_] = s2mpj_ii('CS'+str(int(v_['1'])),ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'CS'+str(int(v_['1'])))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X'+str(int(v_['RTN+1']))]])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X'+str(int(v_['2']))]])
        valA = np.append(valA,float(-1.0))
        for I in range(int(v_['2']),int(v_['RTN-1'])+1):
            v_['I+1'] = 1+I
            v_['I-1'] = -1+I
            v_['I+RTN'] = I+v_['RTN']
            [ig,ig_,_] = s2mpj_ii('CS'+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'CS'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(I)]])
            valA = np.append(valA,float(4.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['I+RTN']))]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['I-1']))]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['I+1']))]])
            valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('CS'+str(int(v_['RTN'])),ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'CS'+str(int(v_['RTN'])))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X'+str(int(v_['RTN']))]])
        valA = np.append(valA,float(4.0))
        [ig,ig_,_] = s2mpj_ii('CS'+str(int(v_['RTN'])),ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'CS'+str(int(v_['RTN'])))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X'+str(int(v_['RTN-1']))]])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X'+str(int(v_['2RTN']))]])
        valA = np.append(valA,float(-1.0))
        v_['JS'] = v_['RTN']
        for J in range(int(v_['RTN+1']),int(v_['M-RTN+1'])+1,int(v_['RTN'])):
            v_['J+1'] = 1+J
            v_['JS'] = J+v_['RTN-1']
            v_['JS-1'] = -1+v_['JS']
            v_['J-RTN'] = J-v_['RTN']
            v_['J+RTN'] = J+v_['RTN']
            [ig,ig_,_] = s2mpj_ii('CS'+str(J),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'CS'+str(J))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(J)]])
            valA = np.append(valA,float(4.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['J+1']))]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['J-RTN']))]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['J+RTN']))]])
            valA = np.append(valA,float(-1.0))
            for I in range(int(v_['J+1']),int(v_['JS-1'])+1):
                v_['I+1'] = 1+I
                v_['I-1'] = -1+I
                v_['I+RTN'] = I+v_['RTN']
                v_['I-RTN'] = I-v_['RTN']
                [ig,ig_,_] = s2mpj_ii('CS'+str(I),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'CS'+str(I))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(I)]])
                valA = np.append(valA,float(4.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(int(v_['I-1']))]])
                valA = np.append(valA,float(-1.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(int(v_['I+1']))]])
                valA = np.append(valA,float(-1.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(int(v_['I-RTN']))]])
                valA = np.append(valA,float(-1.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(int(v_['I+RTN']))]])
                valA = np.append(valA,float(-1.0))
            v_['JS+RTN'] = v_['JS']+v_['RTN']
            v_['JS-RTN'] = v_['JS']-v_['RTN']
            [ig,ig_,_] = s2mpj_ii('CS'+str(int(v_['JS'])),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'CS'+str(int(v_['JS'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['JS']))]])
            valA = np.append(valA,float(4.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['JS-1']))]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('CS'+str(int(v_['JS'])),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'CS'+str(int(v_['JS'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['JS-RTN']))]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['JS+RTN']))]])
            valA = np.append(valA,float(-1.0))
        v_['K'] = 1+v_['JS']
        for I in range(int(v_['K']),int(v_['M'])+1,int(v_['M'])):
            v_['K+1'] = 1+v_['K']
            v_['K+RTN'] = v_['K']+v_['RTN']
            v_['K-RTN'] = v_['K']-v_['RTN']
            [ig,ig_,_] = s2mpj_ii('CS'+str(int(v_['K'])),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'CS'+str(int(v_['K'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['K']))]])
            valA = np.append(valA,float(4.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['K+1']))]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('CS'+str(int(v_['K'])),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'CS'+str(int(v_['K'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['K-RTN']))]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['K+RTN']))]])
            valA = np.append(valA,float(-1.0))
        v_['K'] = 1+v_['K']
        for I in range(int(v_['K']),int(v_['M'])+1):
            v_['I+1'] = 1+I
            v_['I-1'] = -1+I
            v_['I+RTN'] = I+v_['RTN']
            v_['I-RTN'] = I-v_['RTN']
            [ig,ig_,_] = s2mpj_ii('CS'+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'CS'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(I)]])
            valA = np.append(valA,float(4.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['I-1']))]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['I+1']))]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['I-RTN']))]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['I+RTN']))]])
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
        self.gconst = arrset(self.gconst,ig_['CS'+str(int(v_['1']))],float(0.5))
        self.gconst = arrset(self.gconst,ig_['CS'+str(int(v_['RTN']))],float(0.5))
        v_['K'] = v_['RTN+1']
        for J in range(int(v_['RTN+1']),int(v_['M-RTN+1'])+1,int(v_['RTN'])):
            v_['K'] = 1+v_['K']
            for I in range(int(v_['1']),int(v_['RTN-2'])+1):
                v_['K'] = 1+v_['K']
                self.gconst = arrset(self.gconst,ig_['CS'+str(int(v_['K']))],float(-0.5))
            v_['K'] = 1+v_['K']
        v_['K'] = 1+v_['K']
        for J in range(int(v_['K']),int(v_['M'])+1):
            self.gconst = arrset(self.gconst,ig_['CS'+str(J)],float(-0.5))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(0.5))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'X')
        [it,iet_,_] = s2mpj_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['1']),int(v_['N'])+1):
            ename = 'XSQ'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSQ')
            ielftype = arrset(ielftype,ie,iet_["eSQ"])
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        for I in range(int(v_['1']),int(v_['NNZ'])+1):
            v_['RKI'] = v_['K'+str(I)]
            v_['KI'] = int(np.fix(v_['RKI']))
            v_['I-1'] = -1+I
            for J in range(int(v_['1']),int(v_['I-1'])+1):
                v_['RKJ'] = v_['K'+str(J)]
                v_['KJ'] = int(np.fix(v_['RKJ']))
                ename = 'P'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'en2PR')
                ielftype = arrset(ielftype,ie,iet_["en2PR"])
                vname = 'X'+str(int(v_['KI']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
                posev = np.where(elftv[ielftype[ie]]=='X')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'X'+str(int(v_['KJ']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
                posev = np.where(elftv[ielftype[ie]]=='Y')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['TMP'] = 0.5*v_['D'+str(I)]
            ig = ig_['OBJ']
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['XSQ'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['TMP']))
        for I in range(int(v_['1']),int(v_['NNZ'])+1):
            v_['RKI'] = v_['K'+str(I)]
            v_['KI'] = int(np.fix(v_['RKI']))
            v_['I-1'] = -1+I
            for J in range(int(v_['1']),int(v_['I-1'])+1):
                v_['TMP'] = v_['DY'+str(I)]*v_['Y'+str(J)]
                v_['WIJ'] = v_['TMP']*v_['-2/YN2']
                v_['TMP'] = v_['DY'+str(J)]*v_['Y'+str(I)]
                v_['TMP'] = v_['TMP']*v_['-2/YN2']
                v_['WIJ'] = v_['WIJ']+v_['TMP']
                v_['TMP'] = v_['Y'+str(I)]*v_['Y'+str(J)]
                v_['TMP'] = v_['TMP']*v_['DD']
                v_['WIJ'] = v_['WIJ']+v_['TMP']
                ig = ig_['OBJ']
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['P'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(v_['WIJ']))
            v_['TMP'] = v_['DY'+str(I)]*v_['Y'+str(I)]
            v_['WII'] = v_['TMP']*v_['-2/YN2']
            v_['TMP'] = v_['Y'+str(I)]*v_['Y'+str(I)]
            v_['TMP'] = v_['TMP']*v_['DD/2']
            v_['WII'] = v_['WII']+v_['TMP']
            ig = ig_['OBJ']
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['XSQ'+str(int(v_['KI']))])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['WII']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN(  36, 10,2)   -35.69811798
# LO SOLTN( 900, 30,1)   -509.8245900
# LO SOLTN( 900, 30,2)   -950.8404853
# LO SOLTN( 900, 30,3)   -1896.596722
# LO SOLTN( 900, 60,1)   -504.3600140
# LO SOLTN( 900, 60,2)   -945.1134463
# LO SOLTN( 900, 60,3)   -1890.602184
# LO SOLTN( 900, 90,1)   -498.9518964
# LO SOLTN( 900, 90,2)   -939.2704526
# LO SOLTN( 900, 90,3)   -1884.291256
# LO SOLTN( 900,120,1)   -493.5058050
# LO SOLTN( 900,120,2)   -933.1963138
# LO SOLTN( 900,120,3)   -1877.513644
# LO SOLTN( 900,300,1)   -457.1185630
# LO SOLTN( 900,300,2)   -887.3869230
# LO SOLTN( 900,300,3)   -1819.655008
# LO SOLTN( 900,600,1)   -377.5813314
# LO SOLTN( 900,600,2)   -755.0919955
# LO SOLTN( 900,600,3)   -1597.482277
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle+self.neq,self.m)] = np.zeros((self.nge,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CQLR2-AN-V-V"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSQ(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0,0]*EV_[0,0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[0,0]+EV_[0,0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def en2PR(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0,0]*EV_[1,0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1,0]
            g_[1] = EV_[0,0]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = 1.0
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

