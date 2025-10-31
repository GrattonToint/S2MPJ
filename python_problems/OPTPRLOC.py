from s2mpjlib import *
class  OPTPRLOC(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : OPTPRLOC
#    *********
# 
#    Optimal positioning of a new product in a multiattribute space.
#    Consider a market of M existing products, a set of N consumers
#    in a multiattribute (dim K) space.
# 
#    Source: Test problem 4 in M. Duran & I.E. Grossmann,
#    "An outer approximation algorithm for a class of mixed integer nonlinear
#     programs", Mathematical Programming 36, pp. 307-339, 1986.
# 
#    SIF input: S. Leyffer, October 1997
# 
#    classification = "C-CQQR2-AN-30-30"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'OPTPRLOC'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['1'] = 1
        v_['K'] = 5
        v_['M'] = 10
        v_['N'] = 25
        v_['H'] = 1000.0
        v_['Z1,1'] = 2.26
        v_['Z1,2'] = 5.15
        v_['Z1,3'] = 4.03
        v_['Z1,4'] = 1.74
        v_['Z1,5'] = 4.74
        v_['Z2,1'] = 5.51
        v_['Z2,2'] = 9.01
        v_['Z2,3'] = 3.84
        v_['Z2,4'] = 1.47
        v_['Z2,5'] = 9.92
        v_['Z3,1'] = 4.06
        v_['Z3,2'] = 1.80
        v_['Z3,3'] = 0.71
        v_['Z3,4'] = 9.09
        v_['Z3,5'] = 8.13
        v_['Z4,1'] = 6.30
        v_['Z4,2'] = 0.11
        v_['Z4,3'] = 4.08
        v_['Z4,4'] = 7.29
        v_['Z4,5'] = 4.24
        v_['Z5,1'] = 2.81
        v_['Z5,2'] = 1.65
        v_['Z5,3'] = 8.08
        v_['Z5,4'] = 3.99
        v_['Z5,5'] = 3.51
        v_['Z6,1'] = 4.29
        v_['Z6,2'] = 9.49
        v_['Z6,3'] = 2.24
        v_['Z6,4'] = 9.78
        v_['Z6,5'] = 1.52
        v_['Z7,1'] = 9.76
        v_['Z7,2'] = 3.64
        v_['Z7,3'] = 6.62
        v_['Z7,4'] = 3.66
        v_['Z7,5'] = 9.08
        v_['Z8,1'] = 1.37
        v_['Z8,2'] = 6.99
        v_['Z8,3'] = 7.19
        v_['Z8,4'] = 3.03
        v_['Z8,5'] = 3.39
        v_['Z9,1'] = 8.89
        v_['Z9,2'] = 8.29
        v_['Z9,3'] = 6.05
        v_['Z9,4'] = 7.48
        v_['Z9,5'] = 4.09
        v_['Z10,1'] = 7.42
        v_['Z10,2'] = 4.60
        v_['Z10,3'] = 0.30
        v_['Z10,4'] = 0.97
        v_['Z10,5'] = 8.77
        v_['Z11,1'] = 1.54
        v_['Z11,2'] = 7.06
        v_['Z11,3'] = 0.01
        v_['Z11,4'] = 1.23
        v_['Z11,5'] = 3.11
        v_['Z12,1'] = 7.74
        v_['Z12,2'] = 4.40
        v_['Z12,3'] = 7.93
        v_['Z12,4'] = 5.95
        v_['Z12,5'] = 4.88
        v_['Z13,1'] = 9.94
        v_['Z13,2'] = 5.21
        v_['Z13,3'] = 8.58
        v_['Z13,4'] = 0.13
        v_['Z13,5'] = 4.57
        v_['Z14,1'] = 9.54
        v_['Z14,2'] = 1.57
        v_['Z14,3'] = 9.66
        v_['Z14,4'] = 5.24
        v_['Z14,5'] = 7.90
        v_['Z15,1'] = 7.46
        v_['Z15,2'] = 8.81
        v_['Z15,3'] = 1.67
        v_['Z15,4'] = 6.47
        v_['Z15,5'] = 1.81
        v_['Z16,1'] = 0.56
        v_['Z16,2'] = 8.10
        v_['Z16,3'] = 0.19
        v_['Z16,4'] = 6.11
        v_['Z16,5'] = 6.40
        v_['Z17,1'] = 3.86
        v_['Z17,2'] = 6.68
        v_['Z17,3'] = 6.42
        v_['Z17,4'] = 7.29
        v_['Z17,5'] = 4.66
        v_['Z18,1'] = 2.98
        v_['Z18,2'] = 2.98
        v_['Z18,3'] = 3.03
        v_['Z18,4'] = 0.02
        v_['Z18,5'] = 0.67
        v_['Z19,1'] = 3.61
        v_['Z19,2'] = 7.62
        v_['Z19,3'] = 1.79
        v_['Z19,4'] = 7.80
        v_['Z19,5'] = 9.81
        v_['Z20,1'] = 5.68
        v_['Z20,2'] = 4.24
        v_['Z20,3'] = 4.17
        v_['Z20,4'] = 6.75
        v_['Z20,5'] = 1.08
        v_['Z21,1'] = 5.48
        v_['Z21,2'] = 3.74
        v_['Z21,3'] = 3.34
        v_['Z21,4'] = 6.22
        v_['Z21,5'] = 7.94
        v_['Z22,1'] = 8.13
        v_['Z22,2'] = 8.72
        v_['Z22,3'] = 3.93
        v_['Z22,4'] = 8.80
        v_['Z22,5'] = 8.56
        v_['Z23,1'] = 1.37
        v_['Z23,2'] = 0.54
        v_['Z23,3'] = 1.55
        v_['Z23,4'] = 5.56
        v_['Z23,5'] = 5.85
        v_['Z24,1'] = 8.79
        v_['Z24,2'] = 5.04
        v_['Z24,3'] = 4.83
        v_['Z24,4'] = 6.94
        v_['Z24,5'] = 0.38
        v_['Z25,1'] = 2.66
        v_['Z25,2'] = 4.19
        v_['Z25,3'] = 6.49
        v_['Z25,4'] = 8.04
        v_['Z25,5'] = 1.66
        v_['W1,1'] = 9.57
        v_['W1,2'] = 2.74
        v_['W1,3'] = 9.75
        v_['W1,4'] = 3.96
        v_['W1,5'] = 8.67
        v_['W2,1'] = 8.38
        v_['W2,2'] = 3.93
        v_['W2,3'] = 5.18
        v_['W2,4'] = 5.20
        v_['W2,5'] = 7.82
        v_['W3,1'] = 9.81
        v_['W3,2'] = 0.04
        v_['W3,3'] = 4.21
        v_['W3,4'] = 7.38
        v_['W3,5'] = 4.11
        v_['W4,1'] = 7.41
        v_['W4,2'] = 6.08
        v_['W4,3'] = 5.46
        v_['W4,4'] = 4.86
        v_['W4,5'] = 1.48
        v_['W5,1'] = 9.96
        v_['W5,2'] = 9.13
        v_['W5,3'] = 2.95
        v_['W5,4'] = 8.25
        v_['W5,5'] = 3.58
        v_['W6,1'] = 9.39
        v_['W6,2'] = 4.27
        v_['W6,3'] = 5.09
        v_['W6,4'] = 1.81
        v_['W6,5'] = 7.58
        v_['W7,1'] = 1.88
        v_['W7,2'] = 7.20
        v_['W7,3'] = 6.65
        v_['W7,4'] = 1.74
        v_['W7,5'] = 2.86
        v_['W8,1'] = 4.01
        v_['W8,2'] = 2.67
        v_['W8,3'] = 4.86
        v_['W8,4'] = 2.55
        v_['W8,5'] = 6.91
        v_['W9,1'] = 4.18
        v_['W9,2'] = 1.92
        v_['W9,3'] = 2.60
        v_['W9,4'] = 7.15
        v_['W9,5'] = 2.86
        v_['W10,1'] = 7.81
        v_['W10,2'] = 2.14
        v_['W10,3'] = 9.63
        v_['W10,4'] = 7.61
        v_['W10,5'] = 9.17
        v_['W11,1'] = 8.96
        v_['W11,2'] = 3.47
        v_['W11,3'] = 5.49
        v_['W11,4'] = 4.73
        v_['W11,5'] = 9.43
        v_['W12,1'] = 9.94
        v_['W12,2'] = 1.63
        v_['W12,3'] = 1.23
        v_['W12,4'] = 4.33
        v_['W12,5'] = 7.08
        v_['W13,1'] = 0.31
        v_['W13,2'] = 5.00
        v_['W13,3'] = 0.16
        v_['W13,4'] = 2.52
        v_['W13,5'] = 3.08
        v_['W14,1'] = 6.02
        v_['W14,2'] = 0.92
        v_['W14,3'] = 7.47
        v_['W14,4'] = 9.74
        v_['W14,5'] = 1.76
        v_['W15,1'] = 5.06
        v_['W15,2'] = 4.52
        v_['W15,3'] = 1.89
        v_['W15,4'] = 1.22
        v_['W15,5'] = 9.05
        v_['W16,1'] = 5.92
        v_['W16,2'] = 2.56
        v_['W16,3'] = 7.74
        v_['W16,4'] = 6.96
        v_['W16,5'] = 5.18
        v_['W17,1'] = 6.45
        v_['W17,2'] = 1.52
        v_['W17,3'] = 0.06
        v_['W17,4'] = 5.34
        v_['W17,5'] = 8.47
        v_['W18,1'] = 1.04
        v_['W18,2'] = 1.36
        v_['W18,3'] = 5.99
        v_['W18,4'] = 8.10
        v_['W18,5'] = 5.22
        v_['W19,1'] = 1.40
        v_['W19,2'] = 1.35
        v_['W19,3'] = 0.59
        v_['W19,4'] = 8.58
        v_['W19,5'] = 1.21
        v_['W20,1'] = 6.68
        v_['W20,2'] = 9.48
        v_['W20,3'] = 1.60
        v_['W20,4'] = 6.74
        v_['W20,5'] = 8.92
        v_['W21,1'] = 1.95
        v_['W21,2'] = 0.46
        v_['W21,3'] = 2.90
        v_['W21,4'] = 1.79
        v_['W21,5'] = 0.99
        v_['W22,1'] = 5.18
        v_['W22,2'] = 5.10
        v_['W22,3'] = 8.81
        v_['W22,4'] = 3.27
        v_['W22,5'] = 9.63
        v_['W23,1'] = 1.47
        v_['W23,2'] = 5.71
        v_['W23,3'] = 6.95
        v_['W23,4'] = 1.42
        v_['W23,5'] = 3.49
        v_['W24,1'] = 5.40
        v_['W24,2'] = 3.12
        v_['W24,3'] = 5.37
        v_['W24,4'] = 6.10
        v_['W24,5'] = 3.71
        v_['W25,1'] = 6.32
        v_['W25,2'] = 0.81
        v_['W25,3'] = 6.12
        v_['W25,4'] = 6.73
        v_['W25,5'] = 7.93
        v_['DEL1,1'] = 0.62
        v_['DEL1,2'] = 5.06
        v_['DEL1,3'] = 7.82
        v_['DEL1,4'] = 0.22
        v_['DEL1,5'] = 4.42
        v_['DEL2,1'] = 5.21
        v_['DEL2,2'] = 2.66
        v_['DEL2,3'] = 9.54
        v_['DEL2,4'] = 5.03
        v_['DEL2,5'] = 8.01
        v_['DEL3,1'] = 5.27
        v_['DEL3,2'] = 7.72
        v_['DEL3,3'] = 7.97
        v_['DEL3,4'] = 3.31
        v_['DEL3,5'] = 6.56
        v_['DEL4,1'] = 1.02
        v_['DEL4,2'] = 8.89
        v_['DEL4,3'] = 8.77
        v_['DEL4,4'] = 3.10
        v_['DEL4,5'] = 6.66
        v_['DEL5,1'] = 1.26
        v_['DEL5,2'] = 6.80
        v_['DEL5,3'] = 2.30
        v_['DEL5,4'] = 1.75
        v_['DEL5,5'] = 6.65
        v_['DEL6,1'] = 3.74
        v_['DEL6,2'] = 9.06
        v_['DEL6,3'] = 9.80
        v_['DEL6,4'] = 3.01
        v_['DEL6,5'] = 9.52
        v_['DEL7,1'] = 4.64
        v_['DEL7,2'] = 7.99
        v_['DEL7,3'] = 6.69
        v_['DEL7,4'] = 5.88
        v_['DEL7,5'] = 8.23
        v_['DEL8,1'] = 8.35
        v_['DEL8,2'] = 3.79
        v_['DEL8,3'] = 1.19
        v_['DEL8,4'] = 1.96
        v_['DEL8,5'] = 5.88
        v_['DEL9,1'] = 6.44
        v_['DEL9,2'] = 0.17
        v_['DEL9,3'] = 9.93
        v_['DEL9,4'] = 6.80
        v_['DEL9,5'] = 9.75
        v_['DEL10,1'] = 6.49
        v_['DEL10,2'] = 1.92
        v_['DEL10,3'] = 0.05
        v_['DEL10,4'] = 4.89
        v_['DEL10,5'] = 6.43
        v_['R1'] = 77.83985
        v_['R2'] = 175.9710
        v_['R3'] = 201.8226
        v_['R4'] = 143.9533
        v_['R5'] = 154.3895
        v_['R6'] = 433.3177
        v_['R7'] = 109.0764
        v_['R8'] = 41.59592
        v_['R9'] = 144.0623
        v_['R10'] = 99.83416
        v_['R11'] = 149.1791
        v_['R12'] = 123.8074
        v_['R13'] = 27.22197
        v_['R14'] = 89.92683
        v_['R15'] = 293.0766
        v_['R16'] = 174.3170
        v_['R17'] = 125.1028
        v_['R18'] = 222.8417
        v_['R19'] = 50.48593
        v_['R20'] = 361.1973
        v_['R21'] = 40.32642
        v_['R22'] = 161.8518
        v_['R23'] = 66.85827
        v_['R24'] = 340.5807
        v_['R25'] = 407.5200
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['RPH'+str(I)] = v_['R'+str(I)]+v_['H']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['K'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I))
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('Y'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'Y'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('PROFIT',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Y1']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Y2']])
        valA = np.append(valA,float(-0.2))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Y3']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Y4']])
        valA = np.append(valA,float(-0.2))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Y5']])
        valA = np.append(valA,float(-0.9))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Y6']])
        valA = np.append(valA,float(-0.9))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Y7']])
        valA = np.append(valA,float(-0.1))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Y8']])
        valA = np.append(valA,float(-0.8))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Y9']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Y10']])
        valA = np.append(valA,float(-0.4))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Y11']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Y12']])
        valA = np.append(valA,float(-0.3))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Y13']])
        valA = np.append(valA,float(-0.1))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Y14']])
        valA = np.append(valA,float(-0.3))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Y15']])
        valA = np.append(valA,float(-0.5))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Y16']])
        valA = np.append(valA,float(-0.9))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Y17']])
        valA = np.append(valA,float(-0.8))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Y18']])
        valA = np.append(valA,float(-0.1))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Y19']])
        valA = np.append(valA,float(-0.9))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Y20']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Y21']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Y22']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Y23']])
        valA = np.append(valA,float(-0.2))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Y24']])
        valA = np.append(valA,float(-0.7))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Y25']])
        valA = np.append(valA,float(-0.7))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2']])
        valA = np.append(valA,float(-0.9))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X3']])
        valA = np.append(valA,float(-0.5))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X5']])
        valA = np.append(valA,float(1.0))
        for I in range(int(v_['1']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('ELLI'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'ELLI'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['Y'+str(I)]])
            valA = np.append(valA,float(v_['H']))
        [ig,ig_,_] = s2mpj_ii('LIN1',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'LIN1')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X1']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X3']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X4']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X5']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('LIN2',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'LIN2')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X1']])
        valA = np.append(valA,float(0.6))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2']])
        valA = np.append(valA,float(-0.9))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X3']])
        valA = np.append(valA,float(-0.5))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X4']])
        valA = np.append(valA,float(0.1))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X5']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('LIN3',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'LIN3')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X1']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X3']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X4']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X5']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('LIN4',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'LIN4')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X1']])
        valA = np.append(valA,float(0.157))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2']])
        valA = np.append(valA,float(0.05))
        [ig,ig_,_] = s2mpj_ii('LIN5',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'LIN5')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2']])
        valA = np.append(valA,float(0.25))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X4']])
        valA = np.append(valA,float(1.05))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X5']])
        valA = np.append(valA,float(-0.3))
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
        for I in range(int(v_['1']),int(v_['N'])+1):
            self.gconst = arrset(self.gconst,ig_['ELLI'+str(I)],float(v_['RPH'+str(I)]))
        self.gconst = arrset(self.gconst,ig_['LIN1'],float(10.0))
        self.gconst = arrset(self.gconst,ig_['LIN2'],float(-0.64))
        self.gconst = arrset(self.gconst,ig_['LIN3'],float(0.69))
        self.gconst = arrset(self.gconst,ig_['LIN4'],float(1.5))
        self.gconst = arrset(self.gconst,ig_['LIN5'],float(4.5))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        self.xlower[ix_['X1']] = 2.0
        self.xupper[ix_['X1']] = 4.5
        self.xupper[ix_['X2']] = 8.0
        self.xlower[ix_['X3']] = 3.0
        self.xupper[ix_['X3']] = 9.0
        self.xupper[ix_['X4']] = 5.0
        self.xlower[ix_['X5']] = 4.0
        self.xupper[ix_['X5']] = 10.0
        for I in range(int(v_['1']),int(v_['N'])+1):
            self.xupper[ix_['Y'+str(I)]] = 1.0
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eXMBS', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftp = []
        elftp = loaset(elftp,it,0,'B')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        ename = 'OBJSQ1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eXMBS')
        ielftype = arrset(ielftype,ie,iet_["eXMBS"])
        self.x0 = np.zeros((self.n,1))
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='B')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(0.0))
        ename = 'OBJSQ2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eXMBS')
        ielftype = arrset(ielftype,ie,iet_["eXMBS"])
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='B')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(0.0))
        for I in range(int(v_['1']),int(v_['N'])+1):
            for L in range(int(v_['1']),int(v_['K'])+1):
                ename = 'SQ'+str(I)+','+str(L)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eXMBS')
                ielftype = arrset(ielftype,ie,iet_["eXMBS"])
                vname = 'X'+str(L)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                posev = np.where(elftv[ielftype[ie]]=='X')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                posep = np.where(elftp[ielftype[ie]]=='B')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['Z'+str(I)+','+str(L)]))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        ig = ig_['PROFIT']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['OBJSQ1'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.6))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['OBJSQ2'])
        self.grelw = loaset(self.grelw,ig,posel,float(0.1))
        for I in range(int(v_['1']),int(v_['N'])+1):
            for L in range(int(v_['1']),int(v_['K'])+1):
                ig = ig_['ELLI'+str(I)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['SQ'+str(I)+','+str(L)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(v_['W'+str(I)+','+str(L)]))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.cupper[np.arange(self.nle)] = np.zeros((self.nle,1))
        self.clower[np.arange(self.nle+self.neq,self.m)] = np.zeros((self.nge,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CQQR2-AN-30-30"
        self.x0        = np.zeros((self.n,1))
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eXMBS(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = (EV_[0]-self.elpar[iel_][0])**2
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*(EV_[0]-self.elpar[iel_][0])
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

