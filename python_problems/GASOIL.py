from s2mpjlib import *
class  GASOIL(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : GASOIL
#    *********
# 
#    Determine the reaction coefficients for the catalytic cracking of gas oil
#    and other byproducts. The nonlinear model that describes the process is
# 
#      y_1' = - (theta_1 + theta_3 ) y_1^2
#      y_2' = theta_1 y_1^2 + theta_2 y_2
# 
#    with given initial conditions. The problem is to minimize
# 
#     sum{i=1,20} || y(tau_i,theta) - z_i||^2
# 
#    where the z_i are concentration measurements for y at times tau_i (i=1,20)
# 
#    This is problem 12 in the COPS (Version 2) collection of 
#    E. Dolan and J. More'
#    see "Benchmarking Optimization Software with COPS"
#    Argonne National Labs Technical Report ANL/MCS-246 (2000)
# 
#    SIF input: Nick Gould, November 2000
# 
#    classification = "C-COOR2-AN-V-V"
# 
#  The number of differential equations
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'GASOIL'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['NE'] = 2
        if nargin<1:
            v_['NH'] = int(10);  #  SIF file default value
        else:
            v_['NH'] = int(args[0])
#           Alternative values for the SIF file parameters:
# IE NH                  50             $-PARAMETER
# IE NH                  100            $-PARAMETER
# IE NH                  200            $-PARAMETER
# IE NH                  400            $-PARAMETER
        v_['NP'] = 3
        v_['NM'] = 21
        v_['NC'] = 4
        v_['RHO1'] = 0.0694318442
        v_['RHO2'] = 0.3300094782
        v_['RHO3'] = 0.6699905218
        v_['RHO4'] = 0.9305681558
        v_['TAU1'] = 0.0
        v_['TAU2'] = 0.025
        v_['TAU3'] = 0.05
        v_['TAU4'] = 0.075
        v_['TAU5'] = 0.10
        v_['TAU6'] = 0.125
        v_['TAU7'] = 0.150
        v_['TAU8'] = 0.175
        v_['TAU9'] = 0.20
        v_['TAU10'] = 0.225
        v_['TAU11'] = 0.250
        v_['TAU12'] = 0.30
        v_['TAU13'] = 0.35
        v_['TAU14'] = 0.40
        v_['TAU15'] = 0.45
        v_['TAU16'] = 0.50
        v_['TAU17'] = 0.55
        v_['TAU18'] = 0.65
        v_['TAU19'] = 0.75
        v_['TAU20'] = 0.85
        v_['TAU21'] = 0.95
        v_['TF'] = v_['TAU'+str(int(v_['NM']))]
        v_['RNH'] = float(v_['NH'])
        v_['H'] = v_['TF']/v_['RNH']
        v_['Z1,1'] = 1.0000
        v_['Z1,2'] = 0.0000
        v_['Z2,1'] = 0.8105
        v_['Z2,2'] = 0.2000
        v_['Z3,1'] = 0.6208
        v_['Z3,2'] = 0.2886
        v_['Z4,1'] = 0.5258
        v_['Z4,2'] = 0.3010
        v_['Z5,1'] = 0.4345
        v_['Z5,2'] = 0.3215
        v_['Z6,1'] = 0.3903
        v_['Z6,2'] = 0.3123
        v_['Z7,1'] = 0.3342
        v_['Z7,2'] = 0.2716
        v_['Z8,1'] = 0.3034
        v_['Z8,2'] = 0.2551
        v_['Z9,1'] = 0.2735
        v_['Z9,2'] = 0.2258
        v_['Z10,1'] = 0.2405
        v_['Z10,2'] = 0.1959
        v_['Z11,1'] = 0.2283
        v_['Z11,2'] = 0.1789
        v_['Z12,1'] = 0.2071
        v_['Z12,2'] = 0.1457
        v_['Z13,1'] = 0.1669
        v_['Z13,2'] = 0.1198
        v_['Z14,1'] = 0.1530
        v_['Z14,2'] = 0.0909
        v_['Z15,1'] = 0.1339
        v_['Z15,2'] = 0.0719
        v_['Z16,1'] = 0.1265
        v_['Z16,2'] = 0.0561
        v_['Z17,1'] = 0.1200
        v_['Z17,2'] = 0.0460
        v_['Z18,1'] = 0.0990
        v_['Z18,2'] = 0.0280
        v_['Z19,1'] = 0.0870
        v_['Z19,2'] = 0.0190
        v_['Z20,1'] = 0.0770
        v_['Z20,2'] = 0.0140
        v_['Z21,1'] = 0.0690
        v_['Z21,2'] = 0.0100
        v_['BC1'] = 1.0
        v_['BC2'] = 0.0
        v_['1'] = 1
        v_['2'] = 2
        v_['3'] = 3
        v_['NH-1'] = -1+v_['NH']
        v_['FACT0'] = 1.0
        for I in range(int(v_['1']),int(v_['NC'])+1):
            v_['RI'] = float(I)
            v_['I-1'] = -1+I
            v_['FACT'+str(I)] = v_['FACT'+str(int(v_['I-1']))]*v_['RI']
        for I in range(int(v_['1']),int(v_['NM'])+1):
            v_['TAU/H'] = v_['TAU'+str(I)]/v_['H']
            v_['IT/H'] = int(np.fix(v_['TAU/H']))
            v_['IT/H+1'] = 1+v_['IT/H']
            v_['A'] = v_['IT/H+1']
            v_['B'] = v_['NH']
            v_['A'] = -1*v_['A']
            v_['B'] = -1*v_['B']
            v_['A'] = float(v_['A'])
            v_['ABSA'] = np.absolute(v_['A'])
            v_['ABSA'] = int(np.fix(v_['ABSA']))
            v_['B'] = float(v_['B'])
            v_['ABSB'] = np.absolute(v_['B'])
            v_['ABSB'] = int(np.fix(v_['ABSB']))
            v_['ABSA+B'] = v_['ABSA']+v_['ABSB']
            v_['A'] = v_['A']+v_['ABSA+B']
            v_['B'] = v_['B']+v_['ABSA+B']
            v_['A/B'] = int(np.fix(v_['A']/v_['B']))
            v_['B/A'] = int(np.fix(v_['B']/v_['A']))
            v_['SUM'] = v_['A/B']+v_['B/A']
            v_['A'] = v_['A']*v_['A/B']
            v_['B'] = v_['B']*v_['B/A']
            v_['MAXA,B'] = v_['A']+v_['B']
            v_['MAXA,B'] = int(np.fix(v_['MAXA,B']/v_['SUM']))
            v_['MINA,B'] = v_['ABSA+B']-v_['MAXA,B']
            v_['ITAU'+str(I)] = float(v_['MINA,B'])
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['NP'])+1):
            [iv,ix_,_] = s2mpj_ii('THETA'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'THETA'+str(I))
        for I in range(int(v_['1']),int(v_['NH'])+1):
            for J in range(int(v_['1']),int(v_['NE'])+1):
                [iv,ix_,_] = s2mpj_ii('V'+str(I)+','+str(J),ix_)
                self.xnames=arrset(self.xnames,iv,'V'+str(I)+','+str(J))
            for K in range(int(v_['1']),int(v_['NC'])+1):
                for S in range(int(v_['1']),int(v_['NE'])+1):
                    [iv,ix_,_] = s2mpj_ii('W'+str(I)+','+str(K)+','+str(S),ix_)
                    self.xnames=arrset(self.xnames,iv,'W'+str(I)+','+str(K)+','+str(S))
            for J in range(int(v_['1']),int(v_['NC'])+1):
                for S in range(int(v_['1']),int(v_['NE'])+1):
                    [iv,ix_,_] = s2mpj_ii('U'+str(I)+','+str(J)+','+str(S),ix_)
                    self.xnames=arrset(self.xnames,iv,'U'+str(I)+','+str(J)+','+str(S))
                    [iv,ix_,_] = s2mpj_ii('DU'+str(I)+','+str(J)+','+str(S),ix_)
                    self.xnames=arrset(self.xnames,iv,'DU'+str(I)+','+str(J)+','+str(S))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for J in range(int(v_['1']),int(v_['NM'])+1):
            v_['RITAU'] = v_['ITAU'+str(J)]
            v_['I'] = int(np.fix(v_['RITAU']))
            v_['T'] = -1.0+v_['RITAU']
            v_['T'] = v_['T']*v_['H']
            v_['DIFF'] = v_['TAU'+str(J)]-v_['T']
            for S in range(int(v_['1']),int(v_['NE'])+1):
                v_['RATIO'] = v_['DIFF']
                [ig,ig_,_] = s2mpj_ii('OBJ'+str(J)+','+str(S),ig_)
                gtype = arrset(gtype,ig,'<>')
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['V'+str(int(v_['I']))+','+str(S)]])
                valA = np.append(valA,float(1.0))
                for K in range(int(v_['1']),int(v_['NC'])+1):
                    v_['COEF'] = v_['RATIO']/v_['FACT'+str(K)]
                    [ig,ig_,_] = s2mpj_ii('OBJ'+str(J)+','+str(S),ig_)
                    gtype = arrset(gtype,ig,'<>')
                    irA  = np.append(irA,[ig])
                    icA  = np.append(icA,[ix_['W'+str(int(v_['I']))+','+str(K)+','+str(S)]])
                    valA = np.append(valA,float(v_['COEF']))
                    v_['RATIO'] = v_['RATIO']*v_['DIFF']
                    v_['RATIO'] = v_['RATIO']/v_['H']
        for I in range(int(v_['1']),int(v_['NH'])+1):
            for J in range(int(v_['1']),int(v_['NC'])+1):
                v_['RH'] = v_['RHO'+str(J)]
                [ig,ig_,_] = s2mpj_ii('U'+str(I)+','+str(J)+','+str(int(v_['1'])),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'U'+str(I)+','+str(J)+','+str(int(v_['1'])))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['U'+str(I)+','+str(J)+','+str(int(v_['1']))]])
                valA = np.append(valA,float(-1.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['V'+str(I)+','+str(int(v_['1']))]])
                valA = np.append(valA,float(1.0))
                [ig,ig_,_] = s2mpj_ii('U'+str(I)+','+str(J)+','+str(int(v_['2'])),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'U'+str(I)+','+str(J)+','+str(int(v_['2'])))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['U'+str(I)+','+str(J)+','+str(int(v_['2']))]])
                valA = np.append(valA,float(-1.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['V'+str(I)+','+str(int(v_['2']))]])
                valA = np.append(valA,float(1.0))
                v_['PROD'] = v_['RH']*v_['H']
                for K in range(int(v_['1']),int(v_['NC'])+1):
                    v_['COEF'] = v_['PROD']/v_['FACT'+str(K)]
                    [ig,ig_,_] = s2mpj_ii('U'+str(I)+','+str(J)+','+str(int(v_['1'])),ig_)
                    gtype = arrset(gtype,ig,'==')
                    cnames = arrset(cnames,ig,'U'+str(I)+','+str(J)+','+str(int(v_['1'])))
                    irA  = np.append(irA,[ig])
                    icA  = np.append(icA,[ix_['W'+str(I)+','+str(K)+','+str(int(v_['1']))]])
                    valA = np.append(valA,float(v_['COEF']))
                    [ig,ig_,_] = s2mpj_ii('U'+str(I)+','+str(J)+','+str(int(v_['2'])),ig_)
                    gtype = arrset(gtype,ig,'==')
                    cnames = arrset(cnames,ig,'U'+str(I)+','+str(J)+','+str(int(v_['2'])))
                    irA  = np.append(irA,[ig])
                    icA  = np.append(icA,[ix_['W'+str(I)+','+str(K)+','+str(int(v_['2']))]])
                    valA = np.append(valA,float(v_['COEF']))
                    v_['PROD'] = v_['PROD']*v_['RH']
                [ig,ig_,_] = s2mpj_ii('DU'+str(I)+','+str(J)+','+str(int(v_['1'])),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'DU'+str(I)+','+str(J)+','+str(int(v_['1'])))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['DU'+str(I)+','+str(J)+','+str(int(v_['1']))]])
                valA = np.append(valA,float(-1.0))
                [ig,ig_,_] = s2mpj_ii('DU'+str(I)+','+str(J)+','+str(int(v_['2'])),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'DU'+str(I)+','+str(J)+','+str(int(v_['2'])))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['DU'+str(I)+','+str(J)+','+str(int(v_['2']))]])
                valA = np.append(valA,float(-1.0))
                v_['PROD'] = 1.0
                for K in range(int(v_['1']),int(v_['NC'])+1):
                    v_['K-1'] = -1+K
                    v_['COEF'] = v_['PROD']/v_['FACT'+str(int(v_['K-1']))]
                    [ig,ig_,_] = s2mpj_ii('DU'+str(I)+','+str(J)+','+str(int(v_['1'])),ig_)
                    gtype = arrset(gtype,ig,'==')
                    cnames = arrset(cnames,ig,'DU'+str(I)+','+str(J)+','+str(int(v_['1'])))
                    irA  = np.append(irA,[ig])
                    icA  = np.append(icA,[ix_['W'+str(I)+','+str(K)+','+str(int(v_['1']))]])
                    valA = np.append(valA,float(v_['COEF']))
                    [ig,ig_,_] = s2mpj_ii('DU'+str(I)+','+str(J)+','+str(int(v_['2'])),ig_)
                    gtype = arrset(gtype,ig,'==')
                    cnames = arrset(cnames,ig,'DU'+str(I)+','+str(J)+','+str(int(v_['2'])))
                    irA  = np.append(irA,[ig])
                    icA  = np.append(icA,[ix_['W'+str(I)+','+str(K)+','+str(int(v_['2']))]])
                    valA = np.append(valA,float(v_['COEF']))
                    v_['PROD'] = v_['PROD']*v_['RH']
        for I in range(int(v_['1']),int(v_['NH-1'])+1):
            v_['I+1'] = 1+I
            for S in range(int(v_['1']),int(v_['NE'])+1):
                [ig,ig_,_] = s2mpj_ii('C'+str(I)+','+str(S),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'C'+str(I)+','+str(S))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['V'+str(I)+','+str(S)]])
                valA = np.append(valA,float(1.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['V'+str(int(v_['I+1']))+','+str(S)]])
                valA = np.append(valA,float(-1.0))
                for J in range(int(v_['1']),int(v_['NC'])+1):
                    v_['COEF'] = v_['H']/v_['FACT'+str(J)]
                    [ig,ig_,_] = s2mpj_ii('C'+str(I)+','+str(S),ig_)
                    gtype = arrset(gtype,ig,'==')
                    cnames = arrset(cnames,ig,'C'+str(I)+','+str(S))
                    irA  = np.append(irA,[ig])
                    icA  = np.append(icA,[ix_['W'+str(I)+','+str(J)+','+str(S)]])
                    valA = np.append(valA,float(v_['COEF']))
        for I in range(int(v_['1']),int(v_['NH'])+1):
            for J in range(int(v_['1']),int(v_['NC'])+1):
                for S in range(int(v_['1']),int(v_['NE'])+1):
                    [ig,ig_,_] = s2mpj_ii('CO'+str(I)+','+str(J)+','+str(S),ig_)
                    gtype = arrset(gtype,ig,'==')
                    cnames = arrset(cnames,ig,'CO'+str(I)+','+str(J)+','+str(S))
                    irA  = np.append(irA,[ig])
                    icA  = np.append(icA,[ix_['DU'+str(I)+','+str(J)+','+str(S)]])
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
        for J in range(int(v_['1']),int(v_['NM'])+1):
            for S in range(int(v_['1']),int(v_['NE'])+1):
                self.gconst  = (
                      arrset(self.gconst,ig_['OBJ'+str(J)+','+str(S)],float(v_['Z'+str(J)+','+str(S)])))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        for I in range(int(v_['1']),int(v_['NP'])+1):
            self.xlower[ix_['THETA'+str(I)]] = 0.0
        for S in range(int(v_['1']),int(v_['NE'])+1):
            self.xlower[ix_['V'+str(int(v_['1']))+','+str(S)]] = v_['BC'+str(S)]
            self.xupper[ix_['V'+str(int(v_['1']))+','+str(S)]] = v_['BC'+str(S)]
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        for I in range(int(v_['1']),int(v_['NP'])+1):
            if('THETA'+str(I) in ix_):
                self.x0[ix_['THETA'+str(I)]] = float(0.0)
            else:
                self.y0  = (
                      arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['THETA'+str(I)]),float(0.0)))
        for I in range(int(v_['1']),int(v_['NH'])+1):
            for J in range(int(v_['1']),int(v_['NE'])+1):
                if('V'+str(I)+','+str(J) in ix_):
                    self.x0[ix_['V'+str(I)+','+str(J)]] = float(0.0)
                else:
                    self.y0  = (
                          arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['V'+str(I)+','+str(J)]),float(0.0)))
        v_['I1'] = 1
        v_['RITAU'] = v_['ITAU'+str(int(v_['1']))]
        v_['I2'] = int(np.fix(v_['RITAU']))
        for I in range(int(v_['I1']),int(v_['I2'])+1):
            for S in range(int(v_['1']),int(v_['NE'])+1):
                if('V'+str(I)+','+str(S) in ix_):
                    self.x0[ix_['V'+str(I)+','+str(S)]] = float(v_['BC'+str(S)])
                else:
                    self.y0  = (
                          arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['V'+str(I)+','+str(S)]),float(v_['BC'+str(S)])))
                for J in range(int(v_['1']),int(v_['NC'])+1):
                    if('W'+str(I)+','+str(J)+','+str(S) in ix_):
                        self.x0[ix_['W'+str(I)+','+str(J)+','+str(S)]] = float(0.0)
                    else:
                        self.y0  = (
                              arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['W'+str(I)+','+str(J)+','+str(S)]),float(0.0)))
                    if('U'+str(I)+','+str(J)+','+str(S) in ix_):
                        self.x0[ix_['U'+str(I)+','+str(J)+','+str(S)]] = float(v_['BC'+str(S)])
                    else:
                        self.y0  = (
                              arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['U'+str(I)+','+str(J)+','+str(S)]),float(v_['BC'+str(S)])))
                    if('DU'+str(I)+','+str(J)+','+str(S) in ix_):
                        self.x0[ix_['DU'+str(I)+','+str(J)+','+str(S)]] = float(0.0)
                    else:
                        self.y0  = (
                              arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['DU'+str(I)+','+str(J)+','+str(S)]),float(0.0)))
        for K in range(int(v_['2']),int(v_['NM'])+1):
            v_['I1'] = 1+v_['I2']
            v_['RITAU'] = v_['ITAU'+str(K)]
            v_['I2'] = int(np.fix(v_['RITAU']))
            for I in range(int(v_['I1']),int(v_['I2'])+1):
                v_['S'] = 1
                if('V'+str(I)+','+str(int(v_['S'])) in ix_):
                    self.x0[ix_['V'+str(I)+','+str(int(v_['S']))]]  = (
                          float(v_['Z'+str(K)+','+str(int(v_['S']))]))
                else:
                    self.y0  = (
                          arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['V'+str(I)+','+str(int(v_['S']))]),float(v_['Z'+str(K)+','+str(int(v_['S']))])))
                for J in range(int(v_['1']),int(v_['NC'])+1):
                    if('W'+str(I)+','+str(J)+','+str(int(v_['S'])) in ix_):
                        self.x0[ix_['W'+str(I)+','+str(J)+','+str(int(v_['S']))]] = float(0.0)
                    else:
                        self.y0  = (
                              arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['W'+str(I)+','+str(J)+','+str(int(v_['S']))]),float(0.0)))
                    if('U'+str(I)+','+str(J)+','+str(int(v_['S'])) in ix_):
                        self.x0[ix_['U'+str(I)+','+str(J)+','+str(int(v_['S']))]]  = (
                              float(v_['Z'+str(K)+','+str(int(v_['S']))]))
                    else:
                        self.y0  = (
                              arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['U'+str(I)+','+str(J)+','+str(int(v_['S']))]),float(v_['Z'+str(K)+','+str(int(v_['S']))])))
                    if('DU'+str(I)+','+str(J)+','+str(int(v_['S'])) in ix_):
                        self.x0[ix_['DU'+str(I)+','+str(J)+','+str(int(v_['S']))]] = float(0.0)
                    else:
                        self.y0  = (
                              arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['DU'+str(I)+','+str(J)+','+str(int(v_['S']))]),float(0.0)))
                v_['S'] = 2
                if('V'+str(I)+','+str(int(v_['S'])) in ix_):
                    self.x0[ix_['V'+str(I)+','+str(int(v_['S']))]]  = (
                          float(v_['Z'+str(K)+','+str(int(v_['S']))]))
                else:
                    self.y0  = (
                          arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['V'+str(I)+','+str(int(v_['S']))]),float(v_['Z'+str(K)+','+str(int(v_['S']))])))
                for J in range(int(v_['1']),int(v_['NC'])+1):
                    if('W'+str(I)+','+str(J)+','+str(int(v_['S'])) in ix_):
                        self.x0[ix_['W'+str(I)+','+str(J)+','+str(int(v_['S']))]] = float(0.0)
                    else:
                        self.y0  = (
                              arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['W'+str(I)+','+str(J)+','+str(int(v_['S']))]),float(0.0)))
                    if('U'+str(I)+','+str(J)+','+str(int(v_['S'])) in ix_):
                        self.x0[ix_['U'+str(I)+','+str(J)+','+str(int(v_['S']))]]  = (
                              float(v_['Z'+str(K)+','+str(int(v_['S']))]))
                    else:
                        self.y0  = (
                              arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['U'+str(I)+','+str(J)+','+str(int(v_['S']))]),float(v_['Z'+str(K)+','+str(int(v_['S']))])))
                    if('DU'+str(I)+','+str(J)+','+str(int(v_['S'])) in ix_):
                        self.x0[ix_['DU'+str(I)+','+str(J)+','+str(int(v_['S']))]] = float(0.0)
                    else:
                        self.y0  = (
                              arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['DU'+str(I)+','+str(J)+','+str(int(v_['S']))]),float(0.0)))
        v_['I1'] = 1+v_['I2']
        v_['I2'] = v_['NH']
        for I in range(int(v_['I1']),int(v_['I2'])+1):
            for S in range(int(v_['1']),int(v_['NE'])+1):
                if('V'+str(I)+','+str(S) in ix_):
                    self.x0[ix_['V'+str(I)+','+str(S)]]  = (
                          float(v_['Z'+str(int(v_['NM']))+','+str(S)]))
                else:
                    self.y0  = (
                          arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['V'+str(I)+','+str(S)]),float(v_['Z'+str(int(v_['NM']))+','+str(S)])))
                for J in range(int(v_['1']),int(v_['NC'])+1):
                    if('W'+str(I)+','+str(J)+','+str(S) in ix_):
                        self.x0[ix_['W'+str(I)+','+str(J)+','+str(S)]] = float(0.0)
                    else:
                        self.y0  = (
                              arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['W'+str(I)+','+str(J)+','+str(S)]),float(0.0)))
                    if('U'+str(I)+','+str(J)+','+str(S) in ix_):
                        self.x0[ix_['U'+str(I)+','+str(J)+','+str(S)]]  = (
                              float(v_['Z'+str(int(v_['NM']))+','+str(S)]))
                    else:
                        self.y0  = (
                              arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['U'+str(I)+','+str(J)+','+str(S)]),float(v_['Z'+str(int(v_['NM']))+','+str(S)])))
                    if('DU'+str(I)+','+str(J)+','+str(S) in ix_):
                        self.x0[ix_['DU'+str(I)+','+str(J)+','+str(S)]] = float(0.0)
                    else:
                        self.y0  = (
                              arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['DU'+str(I)+','+str(J)+','+str(S)]),float(0.0)))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'ePROD1', iet_)
        elftv = loaset(elftv,it,0,'THETA1')
        elftv = loaset(elftv,it,1,'THETA3')
        elftv = loaset(elftv,it,2,'U')
        [it,iet_,_] = s2mpj_ii( 'ePROD2', iet_)
        elftv = loaset(elftv,it,0,'THETA')
        elftv = loaset(elftv,it,1,'U')
        [it,iet_,_] = s2mpj_ii( 'ePROD3', iet_)
        elftv = loaset(elftv,it,0,'THETA')
        elftv = loaset(elftv,it,1,'U')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['1']),int(v_['NH'])+1):
            for J in range(int(v_['1']),int(v_['NC'])+1):
                ename = 'P1'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'ePROD1')
                ielftype = arrset(ielftype,ie,iet_["ePROD1"])
                vname = 'THETA'+str(int(v_['1']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='THETA1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'THETA'+str(int(v_['3']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='THETA3')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'U'+str(I)+','+str(J)+','+str(int(v_['1']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='U')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                ename = 'P2'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'ePROD2')
                ielftype = arrset(ielftype,ie,iet_["ePROD2"])
                vname = 'THETA'+str(int(v_['1']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='THETA')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'U'+str(I)+','+str(J)+','+str(int(v_['1']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='U')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                ename = 'P3'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'ePROD3')
                ielftype = arrset(ielftype,ie,iet_["ePROD3"])
                vname = 'THETA'+str(int(v_['2']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='THETA')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'U'+str(I)+','+str(J)+','+str(int(v_['2']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='U')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gL2',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['NH'])+1):
            for J in range(int(v_['1']),int(v_['NC'])+1):
                ig = ig_['CO'+str(I)+','+str(J)+','+str(int(v_['1']))]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['P1'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,1.)
                ig = ig_['CO'+str(I)+','+str(J)+','+str(int(v_['2']))]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['P2'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,1.)
                posel = posel+1
                self.grelt = loaset(self.grelt,ig,posel,ie_['P3'+str(I)+','+str(J)])
                self.grelw = loaset(self.grelw,ig,posel, 1.)
        for J in range(int(v_['1']),int(v_['NM'])+1):
            for S in range(int(v_['1']),int(v_['NE'])+1):
                ig = ig_['OBJ'+str(J)+','+str(S)]
                self.grftype = arrset(self.grftype,ig,'gL2')
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLUTION             5.23664D-03   $ (NH=50)
# LO SOLUTION             5.23659D-03   $ (NH=100)
# LO SOLUTION             5.23659D-03   $ (NH=200)
# LO SOLUTION             5.23659D-03   $ (NH=400)
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
        self.pbclass   = "C-COOR2-AN-V-V"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def ePROD1(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((2,3))
        IV_ = np.zeros(2)
        U_[0,0] = U_[0,0]+1
        U_[0,1] = U_[0,1]+1
        U_[1,2] = U_[1,2]+1
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        f_   = IV_[0]*IV_[1]**2
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = IV_[1]**2
            g_[1] = 2.0*IV_[0]*IV_[1]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = 2.0*IV_[1]
                H_[1,0] = H_[0,1]
                H_[1,1] = 2.0*IV_[0]
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ePROD2(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = -EV_[0]*EV_[1]**2
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -EV_[1]**2
            g_[1] = -2.0*EV_[0]*EV_[1]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = -2.0*EV_[1]
                H_[1,0] = H_[0,1]
                H_[1,1] = -2.0*EV_[0]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ePROD3(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]
            g_[1] = EV_[0]
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

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gL2(self,nargout,*args):

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

