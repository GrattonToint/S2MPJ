from s2mpjlib import *
class  SPANHYD(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : SPANHYD
#    *********
# 
#    The spanish hydro-electric reservoir management problem
#    The problem is also named "HYD33" in some references.
#    This is a nonlinear network problem
# 
#    Source:
#    J.L. de la Fuente,
#    "Programacion no-lineal: applicationes en analisis, gestion y
#    planificacion de sistemas electricos",
#    Hidroelectrica Espanola, private communication, 1986.
# 
#    SIF input: Ph. Toint, Sept 1990.
# 
#    classification = "C-CONR2-RN-97-33"
# 
#    Number of arcs    = 97
#    Number of nodes
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'SPANHYD'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['NODES'] = 33
        v_['1'] = 1
        v_['P1A'] = 85.459510678
        v_['P2A'] = -6.255270637
        v_['P3A'] = 1.0862222572
        v_['P1B'] = 81.978824852
        v_['P2B'] = -6.021899007
        v_['P3B'] = 8.4278676858
        v_['P1C'] = 83.420024838
        v_['P2C'] = -6.089102950
        v_['P3C'] = 8.9209611766
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['1']),int(v_['NODES'])+1):
            [ig,ig_,_] = s2mpj_ii('N'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'N'+str(I))
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        ngrp   = len(ig_)
        [iv,ix_,_] = s2mpj_ii('X1',ix_)
        self.xnames=arrset(self.xnames,iv,'X1')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N32']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N1']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X2',ix_)
        self.xnames=arrset(self.xnames,iv,'X2')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N32']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N2']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X3',ix_)
        self.xnames=arrset(self.xnames,iv,'X3')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N32']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N3']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X4',ix_)
        self.xnames=arrset(self.xnames,iv,'X4')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N32']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N4']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X5',ix_)
        self.xnames=arrset(self.xnames,iv,'X5')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N32']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N5']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X6',ix_)
        self.xnames=arrset(self.xnames,iv,'X6')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N32']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N6']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X7',ix_)
        self.xnames=arrset(self.xnames,iv,'X7')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N32']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N7']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X8',ix_)
        self.xnames=arrset(self.xnames,iv,'X8')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N32']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N8']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X9',ix_)
        self.xnames=arrset(self.xnames,iv,'X9')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N32']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N9']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X10',ix_)
        self.xnames=arrset(self.xnames,iv,'X10')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N32']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N10']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X11',ix_)
        self.xnames=arrset(self.xnames,iv,'X11')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N1']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N2']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X12',ix_)
        self.xnames=arrset(self.xnames,iv,'X12')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N2']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N3']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X13',ix_)
        self.xnames=arrset(self.xnames,iv,'X13')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N3']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N5']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X14',ix_)
        self.xnames=arrset(self.xnames,iv,'X14')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N4']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N5']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X15',ix_)
        self.xnames=arrset(self.xnames,iv,'X15')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N5']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N9']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X16',ix_)
        self.xnames=arrset(self.xnames,iv,'X16')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N6']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N7']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X17',ix_)
        self.xnames=arrset(self.xnames,iv,'X17')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N7']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N8']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X18',ix_)
        self.xnames=arrset(self.xnames,iv,'X18')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N8']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N9']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X19',ix_)
        self.xnames=arrset(self.xnames,iv,'X19')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N9']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N10']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X20',ix_)
        self.xnames=arrset(self.xnames,iv,'X20')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N10']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N31']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X21',ix_)
        self.xnames=arrset(self.xnames,iv,'X21')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N1']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N2']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X22',ix_)
        self.xnames=arrset(self.xnames,iv,'X22')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N2']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N3']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X23',ix_)
        self.xnames=arrset(self.xnames,iv,'X23')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N3']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N9']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X24',ix_)
        self.xnames=arrset(self.xnames,iv,'X24')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N4']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N9']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X25',ix_)
        self.xnames=arrset(self.xnames,iv,'X25')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N6']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N7']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X26',ix_)
        self.xnames=arrset(self.xnames,iv,'X26')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N7']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N8']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X27',ix_)
        self.xnames=arrset(self.xnames,iv,'X27')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N8']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N9']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X28',ix_)
        self.xnames=arrset(self.xnames,iv,'X28')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N9']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N10']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X29',ix_)
        self.xnames=arrset(self.xnames,iv,'X29')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N10']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N31']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X30',ix_)
        self.xnames=arrset(self.xnames,iv,'X30')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N1']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N11']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X31',ix_)
        self.xnames=arrset(self.xnames,iv,'X31')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N2']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N12']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X32',ix_)
        self.xnames=arrset(self.xnames,iv,'X32')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N3']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N13']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X33',ix_)
        self.xnames=arrset(self.xnames,iv,'X33')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N4']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N14']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X34',ix_)
        self.xnames=arrset(self.xnames,iv,'X34')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N5']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N15']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X35',ix_)
        self.xnames=arrset(self.xnames,iv,'X35')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N6']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N16']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X36',ix_)
        self.xnames=arrset(self.xnames,iv,'X36')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N7']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N17']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X37',ix_)
        self.xnames=arrset(self.xnames,iv,'X37')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N8']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N18']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X38',ix_)
        self.xnames=arrset(self.xnames,iv,'X38')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N9']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N19']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X39',ix_)
        self.xnames=arrset(self.xnames,iv,'X39')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N10']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N20']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X40',ix_)
        self.xnames=arrset(self.xnames,iv,'X40')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N11']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N12']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X41',ix_)
        self.xnames=arrset(self.xnames,iv,'X41')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N12']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N13']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X42',ix_)
        self.xnames=arrset(self.xnames,iv,'X42')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N13']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N15']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X43',ix_)
        self.xnames=arrset(self.xnames,iv,'X43')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N14']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N15']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X44',ix_)
        self.xnames=arrset(self.xnames,iv,'X44')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N15']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N19']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X45',ix_)
        self.xnames=arrset(self.xnames,iv,'X45')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N16']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N17']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X46',ix_)
        self.xnames=arrset(self.xnames,iv,'X46')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N17']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N18']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X47',ix_)
        self.xnames=arrset(self.xnames,iv,'X47')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N18']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N19']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X48',ix_)
        self.xnames=arrset(self.xnames,iv,'X48')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N19']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N20']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X49',ix_)
        self.xnames=arrset(self.xnames,iv,'X49')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N20']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N31']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X50',ix_)
        self.xnames=arrset(self.xnames,iv,'X50')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N11']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N12']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X51',ix_)
        self.xnames=arrset(self.xnames,iv,'X51')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N12']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N13']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X52',ix_)
        self.xnames=arrset(self.xnames,iv,'X52')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N13']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N19']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X53',ix_)
        self.xnames=arrset(self.xnames,iv,'X53')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N14']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N19']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X54',ix_)
        self.xnames=arrset(self.xnames,iv,'X54')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N16']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N17']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X55',ix_)
        self.xnames=arrset(self.xnames,iv,'X55')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N17']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N18']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X56',ix_)
        self.xnames=arrset(self.xnames,iv,'X56')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N18']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N19']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X57',ix_)
        self.xnames=arrset(self.xnames,iv,'X57')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N19']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N20']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X58',ix_)
        self.xnames=arrset(self.xnames,iv,'X58')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N20']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N31']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X59',ix_)
        self.xnames=arrset(self.xnames,iv,'X59')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N11']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N21']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X60',ix_)
        self.xnames=arrset(self.xnames,iv,'X60')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N12']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N22']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X61',ix_)
        self.xnames=arrset(self.xnames,iv,'X61')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N13']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N23']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X62',ix_)
        self.xnames=arrset(self.xnames,iv,'X62')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N14']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N24']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X63',ix_)
        self.xnames=arrset(self.xnames,iv,'X63')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N15']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N25']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X64',ix_)
        self.xnames=arrset(self.xnames,iv,'X64')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N16']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N26']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X65',ix_)
        self.xnames=arrset(self.xnames,iv,'X65')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N17']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N27']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X66',ix_)
        self.xnames=arrset(self.xnames,iv,'X66')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N18']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N28']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X67',ix_)
        self.xnames=arrset(self.xnames,iv,'X67')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N19']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N29']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X68',ix_)
        self.xnames=arrset(self.xnames,iv,'X68')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N20']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N30']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X69',ix_)
        self.xnames=arrset(self.xnames,iv,'X69')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N21']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N22']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X70',ix_)
        self.xnames=arrset(self.xnames,iv,'X70')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N22']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N23']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X71',ix_)
        self.xnames=arrset(self.xnames,iv,'X71')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N23']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N25']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X72',ix_)
        self.xnames=arrset(self.xnames,iv,'X72')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N24']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N25']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X73',ix_)
        self.xnames=arrset(self.xnames,iv,'X73')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N25']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N29']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X74',ix_)
        self.xnames=arrset(self.xnames,iv,'X74')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N26']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N27']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X75',ix_)
        self.xnames=arrset(self.xnames,iv,'X75')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N27']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N28']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X76',ix_)
        self.xnames=arrset(self.xnames,iv,'X76')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N28']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N29']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X77',ix_)
        self.xnames=arrset(self.xnames,iv,'X77')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N29']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N30']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X78',ix_)
        self.xnames=arrset(self.xnames,iv,'X78')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N30']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N31']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X79',ix_)
        self.xnames=arrset(self.xnames,iv,'X79')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N21']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N22']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X80',ix_)
        self.xnames=arrset(self.xnames,iv,'X80')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N22']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N23']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X81',ix_)
        self.xnames=arrset(self.xnames,iv,'X81')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N23']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N29']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X82',ix_)
        self.xnames=arrset(self.xnames,iv,'X82')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N24']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N29']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X83',ix_)
        self.xnames=arrset(self.xnames,iv,'X83')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N26']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N27']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X84',ix_)
        self.xnames=arrset(self.xnames,iv,'X84')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N27']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N28']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X85',ix_)
        self.xnames=arrset(self.xnames,iv,'X85')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N28']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N29']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X86',ix_)
        self.xnames=arrset(self.xnames,iv,'X86')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N29']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N30']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X87',ix_)
        self.xnames=arrset(self.xnames,iv,'X87')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N30']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N31']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X88',ix_)
        self.xnames=arrset(self.xnames,iv,'X88')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N21']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N33']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X89',ix_)
        self.xnames=arrset(self.xnames,iv,'X89')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N22']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N33']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X90',ix_)
        self.xnames=arrset(self.xnames,iv,'X90')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N23']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N33']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X91',ix_)
        self.xnames=arrset(self.xnames,iv,'X91')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N24']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N33']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X92',ix_)
        self.xnames=arrset(self.xnames,iv,'X92')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N25']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N33']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X93',ix_)
        self.xnames=arrset(self.xnames,iv,'X93')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N26']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N33']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X94',ix_)
        self.xnames=arrset(self.xnames,iv,'X94')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N27']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N33']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X95',ix_)
        self.xnames=arrset(self.xnames,iv,'X95')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N28']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N33']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X96',ix_)
        self.xnames=arrset(self.xnames,iv,'X96')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N29']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N33']])
        valA = np.append(valA,float(1.0))
        [iv,ix_,_] = s2mpj_ii('X97',ix_)
        self.xnames=arrset(self.xnames,iv,'X97')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N30']])
        valA = np.append(valA,float(-1.0))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['N33']])
        valA = np.append(valA,float(1.0))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
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
        self.gconst = arrset(self.gconst,ig_['N1'],float(-5.13800e+01))
        self.gconst = arrset(self.gconst,ig_['N2'],float(-1.38400e+01))
        self.gconst = arrset(self.gconst,ig_['N3'],float(-2.58000))
        self.gconst = arrset(self.gconst,ig_['N4'],float(-2.19100e+01))
        self.gconst = arrset(self.gconst,ig_['N6'],float(-1.29700e+01))
        self.gconst = arrset(self.gconst,ig_['N8'],float(-2.89000))
        self.gconst = arrset(self.gconst,ig_['N9'],float(-2.08400e+01))
        self.gconst = arrset(self.gconst,ig_['N10'],float(-1.71400e+01))
        self.gconst = arrset(self.gconst,ig_['N11'],float(-3.20600e+01))
        self.gconst = arrset(self.gconst,ig_['N12'],float(-2.80000e-01))
        self.gconst = arrset(self.gconst,ig_['N13'],float(-4.20000))
        self.gconst = arrset(self.gconst,ig_['N14'],float(-4.83700e+01))
        self.gconst = arrset(self.gconst,ig_['N16'],float(-1.81300e+01))
        self.gconst = arrset(self.gconst,ig_['N18'],float(1.61000))
        self.gconst = arrset(self.gconst,ig_['N19'],float(-2.66000e+01))
        self.gconst = arrset(self.gconst,ig_['N20'],float(-1.87600e+01))
        self.gconst = arrset(self.gconst,ig_['N21'],float(-1.81300e+01))
        self.gconst = arrset(self.gconst,ig_['N24'],float(-1.81300e+01))
        self.gconst = arrset(self.gconst,ig_['N26'],float(-9.10000))
        self.gconst = arrset(self.gconst,ig_['N28'],float(5.81000))
        self.gconst = arrset(self.gconst,ig_['N29'],float(-9.10000))
        self.gconst = arrset(self.gconst,ig_['N30'],float(-6.02000))
        self.gconst = arrset(self.gconst,ig_['N31'],float(6.08350e+02))
        self.gconst = arrset(self.gconst,ig_['N32'],float(-4.62634e+03))
        self.gconst = arrset(self.gconst,ig_['N33'],float(4.36300e+03))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xupper = np.full((self.n,1),3.02400e+03)
        self.xlower = np.zeros((self.n,1))
        self.xlower[ix_['X1']] = 7.70000e+01
        self.xupper[ix_['X1']] = 7.70100e+01
        self.xlower[ix_['X2']] = 1.12452e+03
        self.xupper[ix_['X2']] = 1.12453e+03
        self.xlower[ix_['X3']] = 1.58000e+02
        self.xupper[ix_['X3']] = 1.58010e+02
        self.xlower[ix_['X4']] = 1.60000e+01
        self.xupper[ix_['X4']] = 1.60100e+01
        self.xlower[ix_['X5']] = 0.00000
        self.xupper[ix_['X5']] = 0.00000
        self.xlower[ix_['X6']] = 7.83650e+02
        self.xupper[ix_['X6']] = 7.83660e+02
        self.xlower[ix_['X7']] = 1.10000e+01
        self.xupper[ix_['X7']] = 1.10100e+01
        self.xlower[ix_['X8']] = 4.90000e+01
        self.xupper[ix_['X8']] = 4.90100e+01
        self.xlower[ix_['X9']] = 2.15517e+03
        self.xupper[ix_['X9']] = 2.15518e+03
        self.xlower[ix_['X10']] = 2.52000e+02
        self.xupper[ix_['X10']] = 2.52010e+02
        self.xupper[ix_['X11']] = 3.97840e+02
        self.xupper[ix_['X12']] = 2.22320e+02
        self.xupper[ix_['X13']] = 2.05630e+02
        self.xupper[ix_['X14']] = 2.05630e+02
        self.xupper[ix_['X15']] = 2.05630e+02
        self.xupper[ix_['X16']] = 1.24830e+02
        self.xupper[ix_['X17']] = 1.27010e+02
        self.xupper[ix_['X18']] = 6.10800e+01
        self.xupper[ix_['X19']] = 6.14840e+02
        self.xupper[ix_['X20']] = 7.78080e+02
        self.xupper[ix_['X25']] = 7.25760e+03
        self.xupper[ix_['X26']] = 1.20960e+03
        self.xupper[ix_['X27']] = 9.07200e+02
        self.xupper[ix_['X28']] = 7.25760e+03
        self.xupper[ix_['X29']] = 7.25760e+03
        self.xlower[ix_['X30']] = 7.70000e+01
        self.xupper[ix_['X30']] = 7.70000e+01
        self.xlower[ix_['X31']] = 4.03400e+02
        self.xupper[ix_['X31']] = 1.31200e+03
        self.xlower[ix_['X32']] = 1.58000e+02
        self.xupper[ix_['X32']] = 1.58000e+02
        self.xlower[ix_['X33']] = 1.60000e+01
        self.xupper[ix_['X33']] = 1.60000e+01
        self.xlower[ix_['X34']] = 0.00000
        self.xupper[ix_['X34']] = 0.00000
        self.xlower[ix_['X35']] = 5.02000e+02
        self.xupper[ix_['X35']] = 9.28460e+02
        self.xlower[ix_['X36']] = 1.10000e+01
        self.xupper[ix_['X36']] = 1.10000e+01
        self.xlower[ix_['X37']] = 4.90000e+01
        self.xupper[ix_['X37']] = 4.90000e+01
        self.xlower[ix_['X38']] = 9.15300e+02
        self.xupper[ix_['X38']] = 2.61160e+03
        self.xlower[ix_['X39']] = 2.52000e+02
        self.xupper[ix_['X39']] = 2.52000e+02
        self.xupper[ix_['X40']] = 3.97840e+02
        self.xupper[ix_['X41']] = 2.22320e+02
        self.xupper[ix_['X42']] = 2.05630e+02
        self.xupper[ix_['X43']] = 2.05630e+02
        self.xupper[ix_['X44']] = 2.05630e+02
        self.xupper[ix_['X45']] = 1.24830e+02
        self.xupper[ix_['X46']] = 1.27010e+02
        self.xupper[ix_['X47']] = 6.10800e+01
        self.xupper[ix_['X48']] = 6.14840e+02
        self.xupper[ix_['X49']] = 7.78080e+02
        self.xupper[ix_['X54']] = 7.25760e+03
        self.xupper[ix_['X55']] = 1.20960e+03
        self.xupper[ix_['X56']] = 9.07200e+02
        self.xupper[ix_['X57']] = 7.25760e+03
        self.xupper[ix_['X58']] = 7.25760e+03
        self.xlower[ix_['X59']] = 7.70000e+01
        self.xupper[ix_['X59']] = 7.70000e+01
        self.xlower[ix_['X60']] = 4.03400e+02
        self.xupper[ix_['X60']] = 1.31200e+03
        self.xlower[ix_['X61']] = 1.58000e+02
        self.xupper[ix_['X61']] = 1.58000e+02
        self.xlower[ix_['X62']] = 1.60000e+01
        self.xupper[ix_['X62']] = 1.60000e+01
        self.xlower[ix_['X63']] = 0.00000
        self.xupper[ix_['X63']] = 0.00000
        self.xlower[ix_['X64']] = 5.05640e+02
        self.xupper[ix_['X64']] = 9.28460e+02
        self.xlower[ix_['X65']] = 1.10000e+01
        self.xupper[ix_['X65']] = 1.10000e+01
        self.xlower[ix_['X66']] = 4.90000e+01
        self.xupper[ix_['X66']] = 4.90000e+01
        self.xlower[ix_['X67']] = 9.15300e+02
        self.xupper[ix_['X67']] = 2.61160e+03
        self.xlower[ix_['X68']] = 2.52000e+02
        self.xupper[ix_['X68']] = 2.52000e+02
        self.xupper[ix_['X69']] = 3.97840e+02
        self.xupper[ix_['X70']] = 2.22320e+02
        self.xupper[ix_['X71']] = 2.05630e+02
        self.xupper[ix_['X72']] = 2.05630e+02
        self.xupper[ix_['X73']] = 2.05630e+02
        self.xupper[ix_['X74']] = 1.24830e+02
        self.xupper[ix_['X75']] = 1.27010e+02
        self.xupper[ix_['X76']] = 6.10800e+01
        self.xupper[ix_['X77']] = 6.14840e+02
        self.xupper[ix_['X78']] = 7.78080e+02
        self.xupper[ix_['X83']] = 7.25760e+03
        self.xupper[ix_['X84']] = 1.20960e+03
        self.xupper[ix_['X85']] = 9.07200e+02
        self.xupper[ix_['X86']] = 7.25760e+03
        self.xupper[ix_['X87']] = 7.25760e+03
        self.xlower[ix_['X88']] = 7.70000e+01
        self.xupper[ix_['X88']] = 7.70100e+01
        self.xlower[ix_['X89']] = 1.10000e+03
        self.xupper[ix_['X89']] = 1.10001e+03
        self.xlower[ix_['X90']] = 1.58000e+02
        self.xupper[ix_['X90']] = 1.58010e+02
        self.xlower[ix_['X91']] = 1.60000e+01
        self.xupper[ix_['X91']] = 1.60100e+01
        self.xlower[ix_['X92']] = 0.00000
        self.xupper[ix_['X92']] = 0.00000
        self.xlower[ix_['X93']] = 7.00000e+02
        self.xupper[ix_['X93']] = 7.00010e+02
        self.xlower[ix_['X94']] = 1.10000e+01
        self.xupper[ix_['X94']] = 1.10100e+01
        self.xlower[ix_['X95']] = 4.90000e+01
        self.xupper[ix_['X95']] = 4.90100e+01
        self.xlower[ix_['X96']] = 2.00000e+03
        self.xupper[ix_['X96']] = 2.00001e+03
        self.xlower[ix_['X97']] = 2.52000e+02
        self.xupper[ix_['X97']] = 2.52010e+02
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(0.0))
        if('X1' in ix_):
            self.x0[ix_['X1']] = float(7.70000e+01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X1']),float(7.70000e+01)))
        if('X2' in ix_):
            self.x0[ix_['X2']] = float(1.12452e+03)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X2'])[0],float(1.12452e+03)))
        if('X3' in ix_):
            self.x0[ix_['X3']] = float(1.58000e+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X3']),float(1.58000e+02)))
        if('X4' in ix_):
            self.x0[ix_['X4']] = float(1.60000e+01)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X4'])[0],float(1.60000e+01)))
        if('X6' in ix_):
            self.x0[ix_['X6']] = float(7.83650e+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X6']),float(7.83650e+02)))
        if('X7' in ix_):
            self.x0[ix_['X7']] = float(1.10000e+01)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X7'])[0],float(1.10000e+01)))
        if('X8' in ix_):
            self.x0[ix_['X8']] = float(4.90000e+01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X8']),float(4.90000e+01)))
        if('X9' in ix_):
            self.x0[ix_['X9']] = float(2.15517e+03)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X9'])[0],float(2.15517e+03)))
        if('X10' in ix_):
            self.x0[ix_['X10']] = float(2.52000e+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X10']),float(2.52000e+02)))
        if('X11' in ix_):
            self.x0[ix_['X11']] = float(5.13800e+01)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X11'])[0],float(5.13800e+01)))
        if('X12' in ix_):
            self.x0[ix_['X12']] = float(1.40210e+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X12']),float(1.40210e+02)))
        if('X13' in ix_):
            self.x0[ix_['X13']] = float(1.42790e+02)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X13'])[0],float(1.42790e+02)))
        if('X14' in ix_):
            self.x0[ix_['X14']] = float(2.19100e+01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X14']),float(2.19100e+01)))
        if('X15' in ix_):
            self.x0[ix_['X15']] = float(1.64700e+02)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X15'])[0],float(1.64700e+02)))
        if('X16' in ix_):
            self.x0[ix_['X16']] = float(5.81900e+01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X16']),float(5.81900e+01)))
        if('X17' in ix_):
            self.x0[ix_['X17']] = float(5.81900e+01)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X17'])[0],float(5.81900e+01)))
        if('X18' in ix_):
            self.x0[ix_['X18']] = float(6.10800e+01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X18']),float(6.10800e+01)))
        if('X19' in ix_):
            self.x0[ix_['X19']] = float(5.66430e+02)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X19'])[0],float(5.66430e+02)))
        if('X20' in ix_):
            self.x0[ix_['X20']] = float(5.83570e+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X20']),float(5.83570e+02)))
        if('X30' in ix_):
            self.x0[ix_['X30']] = float(7.70000e+01)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X30'])[0],float(7.70000e+01)))
        if('X31' in ix_):
            self.x0[ix_['X31']] = float(1.04953e+03)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X31']),float(1.04953e+03)))
        if('X32' in ix_):
            self.x0[ix_['X32']] = float(1.58000e+02)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X32'])[0],float(1.58000e+02)))
        if('X33' in ix_):
            self.x0[ix_['X33']] = float(1.60000e+01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X33']),float(1.60000e+01)))
        if('X35' in ix_):
            self.x0[ix_['X35']] = float(7.38430e+02)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X35'])[0],float(7.38430e+02)))
        if('X36' in ix_):
            self.x0[ix_['X36']] = float(1.10000e+01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X36']),float(1.10000e+01)))
        if('X37' in ix_):
            self.x0[ix_['X37']] = float(4.90000e+01)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X37'])[0],float(4.90000e+01)))
        if('X38' in ix_):
            self.x0[ix_['X38']] = float(1.83536e+03)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X38']),float(1.83536e+03)))
        if('X39' in ix_):
            self.x0[ix_['X39']] = float(2.52000e+02)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X39'])[0],float(2.52000e+02)))
        if('X40' in ix_):
            self.x0[ix_['X40']] = float(3.20600e+01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X40']),float(3.20600e+01)))
        if('X42' in ix_):
            self.x0[ix_['X42']] = float(4.20000)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X42'])[0],float(4.20000)))
        if('X43' in ix_):
            self.x0[ix_['X43']] = float(4.83700e+01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X43']),float(4.83700e+01)))
        if('X44' in ix_):
            self.x0[ix_['X44']] = float(5.25700e+01)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X44'])[0],float(5.25700e+01)))
        if('X45' in ix_):
            self.x0[ix_['X45']] = float(5.98500e+01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X45']),float(5.98500e+01)))
        if('X46' in ix_):
            self.x0[ix_['X46']] = float(5.98500e+01)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X46'])[0],float(5.98500e+01)))
        if('X47' in ix_):
            self.x0[ix_['X47']] = float(5.82400e+01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X47']),float(5.82400e+01)))
        if('X49' in ix_):
            self.x0[ix_['X49']] = float(1.87600e+01)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X49'])[0],float(1.87600e+01)))
        if('X59' in ix_):
            self.x0[ix_['X59']] = float(7.70000e+01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X59']),float(7.70000e+01)))
        if('X60' in ix_):
            self.x0[ix_['X60']] = float(1.08187e+03)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X60'])[0],float(1.08187e+03)))
        if('X61' in ix_):
            self.x0[ix_['X61']] = float(1.58000e+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X61']),float(1.58000e+02)))
        if('X62' in ix_):
            self.x0[ix_['X62']] = float(1.60000e+01)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X62'])[0],float(1.60000e+01)))
        if('X64' in ix_):
            self.x0[ix_['X64']] = float(6.96710e+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X64']),float(6.96710e+02)))
        if('X65' in ix_):
            self.x0[ix_['X65']] = float(1.10000e+01)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X65'])[0],float(1.10000e+01)))
        if('X66' in ix_):
            self.x0[ix_['X66']] = float(4.90000e+01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X66']),float(4.90000e+01)))
        if('X67' in ix_):
            self.x0[ix_['X67']] = float(1.97277e+03)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X67'])[0],float(1.97277e+03)))
        if('X68' in ix_):
            self.x0[ix_['X68']] = float(2.52000e+02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X68']),float(2.52000e+02)))
        if('X69' in ix_):
            self.x0[ix_['X69']] = float(1.81300e+01)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X69'])[0],float(1.81300e+01)))
        if('X72' in ix_):
            self.x0[ix_['X72']] = float(1.81300e+01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X72']),float(1.81300e+01)))
        if('X73' in ix_):
            self.x0[ix_['X73']] = float(1.81300e+01)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X73'])[0],float(1.81300e+01)))
        if('X74' in ix_):
            self.x0[ix_['X74']] = float(5.81000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X74']),float(5.81000)))
        if('X75' in ix_):
            self.x0[ix_['X75']] = float(5.81000)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X75'])[0],float(5.81000)))
        if('X78' in ix_):
            self.x0[ix_['X78']] = float(6.02000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X78']),float(6.02000)))
        if('X88' in ix_):
            self.x0[ix_['X88']] = float(7.70000e+01)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X88'])[0],float(7.70000e+01)))
        if('X89' in ix_):
            self.x0[ix_['X89']] = float(1.10000e+03)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X89']),float(1.10000e+03)))
        if('X90' in ix_):
            self.x0[ix_['X90']] = float(1.58000e+02)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X90'])[0],float(1.58000e+02)))
        if('X91' in ix_):
            self.x0[ix_['X91']] = float(1.60000e+01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X91']),float(1.60000e+01)))
        if('X93' in ix_):
            self.x0[ix_['X93']] = float(7.00000e+02)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X93'])[0],float(7.00000e+02)))
        if('X94' in ix_):
            self.x0[ix_['X94']] = float(1.10000e+01)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X94']),float(1.10000e+01)))
        if('X95' in ix_):
            self.x0[ix_['X95']] = float(4.90000e+01)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X95'])[0],float(4.90000e+01)))
        if('X96' in ix_):
            self.x0[ix_['X96']] = float(2.00000e+03)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X96']),float(2.00000e+03)))
        if('X97' in ix_):
            self.x0[ix_['X97']] = float(2.52000e+02)
        else:
            self.y0  = (
                  arrset(self.y0,np.where(self.congrps==ig_['X97'])[0],float(2.52000e+02)))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSPHYD', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftv = loaset(elftv,it,3,'V4')
        elftv = loaset(elftv,it,4,'V6')
        elftv = loaset(elftv,it,5,'V16')
        elftv = loaset(elftv,it,6,'V7')
        elftv = loaset(elftv,it,7,'V8')
        elftv = loaset(elftv,it,8,'V9')
        elftv = loaset(elftv,it,9,'V10')
        elftv = loaset(elftv,it,10,'V11')
        elftv = loaset(elftv,it,11,'V12')
        elftv = loaset(elftv,it,12,'V13')
        elftv = loaset(elftv,it,13,'V14')
        elftv = loaset(elftv,it,14,'V17')
        elftv = loaset(elftv,it,15,'V18')
        elftv = loaset(elftv,it,16,'V19')
        elftv = loaset(elftv,it,17,'V20')
        elftp = []
        elftp = loaset(elftp,it,0,'P1')
        elftp = loaset(elftp,it,1,'P2')
        elftp = loaset(elftp,it,2,'P3')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        ename = 'EA'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSPHYD')
            ielftype = arrset(ielftype,ie,iet_['eSPHYD'])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V4')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V6')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V7')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X8'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V8')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X9'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V9')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X10'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V10')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X11'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V11')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X12'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V12')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X13'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V13')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X14'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V14')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X16'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V16')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X17'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V17')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X18'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V18')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X19'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V19')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X20'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V20')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='P1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(v_['P1A']))
        posep = np.where(elftp[ielftype[ie]]=='P2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(v_['P2A']))
        posep = np.where(elftp[ielftype[ie]]=='P3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(v_['P3A']))
        ename = 'EB'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSPHYD')
            ielftype = arrset(ielftype,ie,iet_['eSPHYD'])
        vname = 'X30'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X31'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X32'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X33'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V4')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X35'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V6')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X36'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V7')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X37'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V8')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X38'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V9')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X39'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V10')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X40'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V11')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X41'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V12')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X42'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V13')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X43'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V14')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X45'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V16')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X46'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V17')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X47'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V18')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X48'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V19')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X49'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V20')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='P1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(v_['P1B']))
        posep = np.where(elftp[ielftype[ie]]=='P2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(v_['P2B']))
        posep = np.where(elftp[ielftype[ie]]=='P3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(v_['P3B']))
        ename = 'EC'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSPHYD')
            ielftype = arrset(ielftype,ie,iet_['eSPHYD'])
        vname = 'X59'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X60'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X61'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X62'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V4')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X64'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V6')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X65'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V7')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X66'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V8')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X67'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V9')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X68'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V10')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X69'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V11')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X70'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V12')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X71'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V13')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X72'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V14')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X74'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V16')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X75'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V17')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X76'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V18')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X77'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V19')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X78'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(3.02400e+03),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V20')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='P1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(v_['P1C']))
        posep = np.where(elftp[ielftype[ie]]=='P2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(v_['P2C']))
        posep = np.where(elftp[ielftype[ie]]=='P3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(v_['P3C']))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        ig = ig_['OBJ']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['EA'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['EB'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['EC'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               239.738001
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
        self.pbclass   = "C-CONR2-RN-97-33"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def e_globs(self):

        import numpy as np
        self.efpar = np.array([]);
        self.efpar = arrset( self.efpar,0,-637.993)
        self.efpar = arrset( self.efpar,1,-15.4452)
        self.efpar = arrset( self.efpar,2,-6.2597)
        self.efpar = arrset( self.efpar,3,-2.8699)
        self.efpar = arrset( self.efpar,4,-8.6773)
        self.efpar = arrset( self.efpar,5,-8.9261)
        self.efpar = arrset( self.efpar,6,-672.26)
        self.efpar = arrset( self.efpar,7,-334.117)
        self.efpar = arrset( self.efpar,8,-205.321)
        self.efpar = arrset( self.efpar,9,2.02)
        self.efpar = arrset( self.efpar,10,2.57)
        self.efpar = arrset( self.efpar,11,2.2146)
        self.efpar = arrset( self.efpar,12,2.5128)
        self.efpar = arrset( self.efpar,13,2.5707)
        self.efpar = arrset( self.efpar,14,2.7601)
        self.efpar = arrset( self.efpar,15,2.5402)
        self.efpar = arrset( self.efpar,16,2.711)
        self.efpar = arrset( self.efpar,17,2.6673)
        self.efpar = arrset( self.efpar,18,344.904)
        self.efpar = arrset( self.efpar,19,283.67)
        self.efpar = arrset( self.efpar,20,188.597)
        self.efpar = arrset( self.efpar,21,212.864)
        self.efpar = arrset( self.efpar,22,353.204)
        self.efpar = arrset( self.efpar,23,316.895)
        self.efpar = arrset( self.efpar,24,295.055)
        self.efpar = arrset( self.efpar,25,163.013)
        self.efpar = arrset( self.efpar,26,107.338)
        self.efpar = arrset( self.efpar,27,0.1265)
        self.efpar = arrset( self.efpar,28,0.031)
        self.efpar = arrset( self.efpar,29,0.558)
        self.efpar = arrset( self.efpar,30,0.7584)
        self.efpar = arrset( self.efpar,31,0.0541)
        self.efpar = arrset( self.efpar,32,0.9038)
        self.efpar = arrset( self.efpar,33,0.1557)
        self.efpar = arrset( self.efpar,34,0.0262)
        self.efpar = arrset( self.efpar,35,-0.3111)
        self.efpar = arrset( self.efpar,36,-2.5225e-4)
        self.efpar = arrset( self.efpar,37,-7.2e-6)
        self.efpar = arrset( self.efpar,38,-1.4052e-3)
        self.efpar = arrset( self.efpar,39,-1.4353e-2)
        self.efpar = arrset( self.efpar,40,-2.01e-5)
        self.efpar = arrset( self.efpar,41,1.4139e-3)
        self.efpar = arrset( self.efpar,42,1.7195e-3)
        self.efpar = arrset( self.efpar,43,-2.4e-6)
        self.efpar = arrset( self.efpar,44,2.372e-4)
        return pbm

    @staticmethod
    def eSPHYD(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        Z1 = self.efpar[18]+self.efpar[27]*EV_[0]+self.efpar[36]*EV_[0]*EV_[0]
        Z2 = self.efpar[19]+self.efpar[28]*EV_[1]+self.efpar[37]*EV_[1]*EV_[1]
        Z3 = self.efpar[20]+self.efpar[29]*EV_[2]+self.efpar[38]*EV_[2]*EV_[2]
        Z4 = self.efpar[21]+self.efpar[30]*EV_[3]+self.efpar[39]*EV_[3]*EV_[3]
        Z6 = self.efpar[22]+self.efpar[31]*EV_[4]+self.efpar[40]*EV_[4]*EV_[4]
        Z7 = self.efpar[23]+self.efpar[32]*EV_[6]+self.efpar[41]*EV_[6]*EV_[6]
        Z8 = self.efpar[24]+self.efpar[33]*EV_[7]+self.efpar[42]*EV_[7]*EV_[7]
        Z9 = self.efpar[25]+self.efpar[34]*EV_[8]+self.efpar[43]*EV_[8]*EV_[8]
        ZT = self.efpar[26]+self.efpar[35]*EV_[9]+self.efpar[44]*EV_[9]*EV_[9]
        DZ1 = self.efpar[27]+2.0*self.efpar[36]*EV_[0]
        DZ2 = self.efpar[28]+2.0*self.efpar[37]*EV_[1]
        DZ3 = self.efpar[29]+2.0*self.efpar[38]*EV_[2]
        DZ4 = self.efpar[30]+2.0*self.efpar[39]*EV_[3]
        DZ6 = self.efpar[31]+2.0*self.efpar[40]*EV_[4]
        DZ7 = self.efpar[32]+2.0*self.efpar[41]*EV_[6]
        DZ8 = self.efpar[33]+2.0*self.efpar[42]*EV_[7]
        DZ9 = self.efpar[34]+2.0*self.efpar[43]*EV_[8]
        DZT = self.efpar[35]+2.0*self.efpar[44]*EV_[9]
        E1 = self.efpar[0]+self.efpar[9]*Z1
        E2 = self.efpar[1]+self.efpar[10]*(Z2-Z3)
        E3 = self.efpar[2]+self.efpar[11]*(Z3-Z9)
        E4 = self.efpar[3]+self.efpar[12]*(Z4-Z9)
        E6 = self.efpar[4]+self.efpar[13]*(Z6-Z7)
        E7 = self.efpar[5]+self.efpar[14]*(Z7-Z8)
        E8 = self.efpar[6]+self.efpar[15]*Z8
        E9 = self.efpar[7]+self.efpar[16]*Z9
        ET = self.efpar[8]+self.efpar[17]*ZT
        DE11 = self.efpar[9]*DZ1
        DE22 = self.efpar[10]*DZ2
        DE23 = -self.efpar[10]*DZ3
        DE33 = self.efpar[11]*DZ3
        DE39 = -self.efpar[11]*DZ9
        DE44 = self.efpar[12]*DZ4
        DE49 = -self.efpar[12]*DZ9
        DE66 = self.efpar[13]*DZ6
        DE67 = -self.efpar[13]*DZ7
        DE77 = self.efpar[14]*DZ7
        DE78 = -self.efpar[14]*DZ8
        DE88 = self.efpar[15]*DZ8
        DE99 = self.efpar[16]*DZ9
        DETT = self.efpar[17]*DZT
        HE111 = 2.0*self.efpar[9]*self.efpar[36]
        HE222 = 2.0*self.efpar[10]*self.efpar[37]
        HE233 = -2.0*self.efpar[10]*self.efpar[38]
        HE333 = 2.0*self.efpar[11]*self.efpar[38]
        HE399 = -2.0*self.efpar[11]*self.efpar[43]
        HE444 = 2.0*self.efpar[12]*self.efpar[39]
        HE499 = -2.0*self.efpar[12]*self.efpar[43]
        HE666 = 2.0*self.efpar[13]*self.efpar[40]
        HE677 = -2.0*self.efpar[13]*self.efpar[41]
        HE777 = 2.0*self.efpar[14]*self.efpar[41]
        HE788 = -2.0*self.efpar[14]*self.efpar[42]
        HE888 = 2.0*self.efpar[15]*self.efpar[42]
        HE999 = 2.0*self.efpar[16]*self.efpar[43]
        HETTT = 2.0*self.efpar[17]*self.efpar[44]
        PS = (EV_[10]*E1+EV_[11]*E2+EV_[12]*E3+EV_[13]*E4+EV_[5]*E6+EV_[14]*E7+
             EV_[15]*E8+EV_[16]*E9+EV_[17]*ET)
        DP1 = EV_[10]*DE11
        DP2 = EV_[11]*DE22
        DP3 = EV_[12]*DE33+EV_[11]*DE23
        DP4 = EV_[13]*DE44
        DP6 = EV_[5]*DE66
        DP7 = EV_[14]*DE77+EV_[5]*DE67
        DP8 = EV_[15]*DE88+EV_[14]*DE78
        DP9 = EV_[16]*DE99+EV_[12]*DE39+EV_[13]*DE49
        DPT = EV_[17]*DETT
        HP11 = EV_[10]*HE111
        HP22 = EV_[11]*HE222
        HP33 = EV_[12]*HE333+EV_[11]*HE233
        HP44 = EV_[13]*HE444
        HP66 = EV_[5]*HE666
        HP77 = EV_[14]*HE777+EV_[5]*HE677
        HP88 = EV_[15]*HE888+EV_[14]*HE788
        HP99 = EV_[16]*HE999+EV_[12]*HE399+EV_[13]*HE499
        HPTT = EV_[17]*HETTT
        HFPS = self.elpar[iel_][2]+self.elpar[iel_][2]
        DFPS = self.elpar[iel_][1]+HFPS*PS
        f_   = self.elpar[iel_][0]+self.elpar[iel_][1]*PS+self.elpar[iel_][2]*PS*PS
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = DP1*DFPS
            g_[1] = DP2*DFPS
            g_[2] = DP3*DFPS
            g_[3] = DP4*DFPS
            g_[4] = DP6*DFPS
            g_[6] = DP7*DFPS
            g_[7] = DP8*DFPS
            g_[8] = DP9*DFPS
            g_[9] = DPT*DFPS
            g_[10] = E1*DFPS
            g_[11] = E2*DFPS
            g_[12] = E3*DFPS
            g_[13] = E4*DFPS
            g_[5] = E6*DFPS
            g_[14] = E7*DFPS
            g_[15] = E8*DFPS
            g_[16] = E9*DFPS
            g_[17] = ET*DFPS
            if nargout>2:
                H_ = np.zeros((18,18))
                H_[0,0] = HFPS*DP1*DP1+DFPS*HP11
                H_[0,1] = HFPS*DP1*DP2
                H_[1,0] = H_[0,1]
                H_[0,2] = HFPS*DP1*DP3
                H_[2,0] = H_[0,2]
                H_[0,3] = HFPS*DP1*DP4
                H_[3,0] = H_[0,3]
                H_[0,4] = HFPS*DP1*DP6
                H_[4,0] = H_[0,4]
                H_[0,6] = HFPS*DP1*DP7
                H_[6,0] = H_[0,6]
                H_[0,7] = HFPS*DP1*DP8
                H_[7,0] = H_[0,7]
                H_[0,8] = HFPS*DP1*DP9
                H_[8,0] = H_[0,8]
                H_[0,9] = HFPS*DP1*DPT
                H_[9,0] = H_[0,9]
                H_[0,10] = HFPS*DP1*E1+DFPS*DE11
                H_[10,0] = H_[0,10]
                H_[0,11] = HFPS*DP1*E2
                H_[11,0] = H_[0,11]
                H_[0,12] = HFPS*DP1*E3
                H_[12,0] = H_[0,12]
                H_[0,13] = HFPS*DP1*E4
                H_[13,0] = H_[0,13]
                H_[0,5] = HFPS*DP1*E6
                H_[5,0] = H_[0,5]
                H_[0,14] = HFPS*DP1*E7
                H_[14,0] = H_[0,14]
                H_[0,15] = HFPS*DP1*E8
                H_[15,0] = H_[0,15]
                H_[0,16] = HFPS*DP1*E9
                H_[16,0] = H_[0,16]
                H_[0,17] = HFPS*DP1*ET
                H_[17,0] = H_[0,17]
                H_[1,1] = HFPS*DP2*DP2+DFPS*HP22
                H_[1,2] = HFPS*DP2*DP3
                H_[2,1] = H_[1,2]
                H_[1,3] = HFPS*DP2*DP4
                H_[3,1] = H_[1,3]
                H_[1,4] = HFPS*DP2*DP6
                H_[4,1] = H_[1,4]
                H_[1,6] = HFPS*DP2*DP7
                H_[6,1] = H_[1,6]
                H_[1,7] = HFPS*DP2*DP8
                H_[7,1] = H_[1,7]
                H_[1,8] = HFPS*DP2*DP9
                H_[8,1] = H_[1,8]
                H_[1,9] = HFPS*DP2*DPT
                H_[9,1] = H_[1,9]
                H_[1,10] = HFPS*DP2*E1
                H_[10,1] = H_[1,10]
                H_[1,11] = HFPS*DP2*E2+DFPS*DE22
                H_[11,1] = H_[1,11]
                H_[1,12] = HFPS*DP2*E3
                H_[12,1] = H_[1,12]
                H_[1,13] = HFPS*DP2*E4
                H_[13,1] = H_[1,13]
                H_[1,5] = HFPS*DP2*E6
                H_[5,1] = H_[1,5]
                H_[1,14] = HFPS*DP2*E7
                H_[14,1] = H_[1,14]
                H_[1,15] = HFPS*DP2*E8
                H_[15,1] = H_[1,15]
                H_[1,16] = HFPS*DP2*E9
                H_[16,1] = H_[1,16]
                H_[1,17] = HFPS*DP2*ET
                H_[17,1] = H_[1,17]
                H_[2,2] = HFPS*DP3*DP3+DFPS*HP33
                H_[2,3] = HFPS*DP3*DP4
                H_[3,2] = H_[2,3]
                H_[2,4] = HFPS*DP3*DP6
                H_[4,2] = H_[2,4]
                H_[2,6] = HFPS*DP3*DP7
                H_[6,2] = H_[2,6]
                H_[2,7] = HFPS*DP3*DP8
                H_[7,2] = H_[2,7]
                H_[2,8] = HFPS*DP3*DP9
                H_[8,2] = H_[2,8]
                H_[2,9] = HFPS*DP3*DPT
                H_[9,2] = H_[2,9]
                H_[2,10] = HFPS*DP3*E1
                H_[10,2] = H_[2,10]
                H_[2,11] = HFPS*DP3*E2+DFPS*DE23
                H_[11,2] = H_[2,11]
                H_[2,12] = HFPS*DP3*E3+DFPS*DE33
                H_[12,2] = H_[2,12]
                H_[2,13] = HFPS*DP3*E4
                H_[13,2] = H_[2,13]
                H_[2,5] = HFPS*DP3*E6
                H_[5,2] = H_[2,5]
                H_[2,14] = HFPS*DP3*E7
                H_[14,2] = H_[2,14]
                H_[2,15] = HFPS*DP3*E8
                H_[15,2] = H_[2,15]
                H_[2,16] = HFPS*DP3*E9
                H_[16,2] = H_[2,16]
                H_[2,17] = HFPS*DP3*ET
                H_[17,2] = H_[2,17]
                H_[3,3] = HFPS*DP4*DP4+DFPS*HP44
                H_[3,4] = HFPS*DP4*DP6
                H_[4,3] = H_[3,4]
                H_[3,6] = HFPS*DP4*DP7
                H_[6,3] = H_[3,6]
                H_[3,7] = HFPS*DP4*DP8
                H_[7,3] = H_[3,7]
                H_[3,8] = HFPS*DP4*DP9
                H_[8,3] = H_[3,8]
                H_[3,9] = HFPS*DP4*DPT
                H_[9,3] = H_[3,9]
                H_[3,10] = HFPS*DP4*E1
                H_[10,3] = H_[3,10]
                H_[3,11] = HFPS*DP4*E2
                H_[11,3] = H_[3,11]
                H_[3,12] = HFPS*DP4*E3
                H_[12,3] = H_[3,12]
                H_[3,13] = HFPS*DP4*E4+DFPS*DE44
                H_[13,3] = H_[3,13]
                H_[3,5] = HFPS*DP4*E6
                H_[5,3] = H_[3,5]
                H_[3,14] = HFPS*DP4*E7
                H_[14,3] = H_[3,14]
                H_[3,15] = HFPS*DP4*E8
                H_[15,3] = H_[3,15]
                H_[3,16] = HFPS*DP4*E9
                H_[16,3] = H_[3,16]
                H_[3,17] = HFPS*DP4*ET
                H_[17,3] = H_[3,17]
                H_[4,4] = HFPS*DP6*DP6+DFPS*HP66
                H_[4,6] = HFPS*DP6*DP7
                H_[6,4] = H_[4,6]
                H_[4,7] = HFPS*DP6*DP8
                H_[7,4] = H_[4,7]
                H_[4,8] = HFPS*DP6*DP9
                H_[8,4] = H_[4,8]
                H_[4,9] = HFPS*DP6*DPT
                H_[9,4] = H_[4,9]
                H_[4,10] = HFPS*DP6*E1
                H_[10,4] = H_[4,10]
                H_[4,11] = HFPS*DP6*E2
                H_[11,4] = H_[4,11]
                H_[4,12] = HFPS*DP6*E3
                H_[12,4] = H_[4,12]
                H_[4,13] = HFPS*DP6*E4
                H_[13,4] = H_[4,13]
                H_[4,5] = HFPS*DP6*E6+DFPS*DE66
                H_[5,4] = H_[4,5]
                H_[4,14] = HFPS*DP6*E7
                H_[14,4] = H_[4,14]
                H_[4,15] = HFPS*DP6*E8
                H_[15,4] = H_[4,15]
                H_[4,16] = HFPS*DP6*E9
                H_[16,4] = H_[4,16]
                H_[4,17] = HFPS*DP6*ET
                H_[17,4] = H_[4,17]
                H_[6,6] = HFPS*DP7*DP7+DFPS*HP77
                H_[6,7] = HFPS*DP7*DP8
                H_[7,6] = H_[6,7]
                H_[6,8] = HFPS*DP7*DP9
                H_[8,6] = H_[6,8]
                H_[6,9] = HFPS*DP7*DPT
                H_[9,6] = H_[6,9]
                H_[6,10] = HFPS*DP7*E1
                H_[10,6] = H_[6,10]
                H_[6,11] = HFPS*DP7*E2
                H_[11,6] = H_[6,11]
                H_[6,12] = HFPS*DP7*E3
                H_[12,6] = H_[6,12]
                H_[6,13] = HFPS*DP7*E4
                H_[13,6] = H_[6,13]
                H_[6,5] = HFPS*DP7*E6+DFPS*DE67
                H_[5,6] = H_[6,5]
                H_[6,14] = HFPS*DP7*E7+DFPS*DE77
                H_[14,6] = H_[6,14]
                H_[6,15] = HFPS*DP7*E8
                H_[15,6] = H_[6,15]
                H_[6,16] = HFPS*DP7*E9
                H_[16,6] = H_[6,16]
                H_[6,17] = HFPS*DP7*ET
                H_[17,6] = H_[6,17]
                H_[7,7] = HFPS*DP8*DP8+DFPS*HP88
                H_[7,8] = HFPS*DP8*DP9
                H_[8,7] = H_[7,8]
                H_[7,9] = HFPS*DP8*DPT
                H_[9,7] = H_[7,9]
                H_[7,10] = HFPS*DP8*E1
                H_[10,7] = H_[7,10]
                H_[7,11] = HFPS*DP8*E2
                H_[11,7] = H_[7,11]
                H_[7,12] = HFPS*DP8*E3
                H_[12,7] = H_[7,12]
                H_[7,13] = HFPS*DP8*E4
                H_[13,7] = H_[7,13]
                H_[7,5] = HFPS*DP8*E6
                H_[5,7] = H_[7,5]
                H_[7,14] = HFPS*DP8*E7+DFPS*DE78
                H_[14,7] = H_[7,14]
                H_[7,15] = HFPS*DP8*E8+DFPS*DE88
                H_[15,7] = H_[7,15]
                H_[7,16] = HFPS*DP8*E9
                H_[16,7] = H_[7,16]
                H_[7,17] = HFPS*DP8*ET
                H_[17,7] = H_[7,17]
                H_[8,8] = HFPS*DP9*DP9+DFPS*HP99
                H_[8,9] = HFPS*DP9*DPT
                H_[9,8] = H_[8,9]
                H_[8,10] = HFPS*DP9*E1
                H_[10,8] = H_[8,10]
                H_[8,11] = HFPS*DP9*E2
                H_[11,8] = H_[8,11]
                H_[8,12] = HFPS*DP9*E3+DFPS*DE39
                H_[12,8] = H_[8,12]
                H_[8,13] = HFPS*DP9*E4+DFPS*DE49
                H_[13,8] = H_[8,13]
                H_[8,5] = HFPS*DP9*E6
                H_[5,8] = H_[8,5]
                H_[8,14] = HFPS*DP9*E7
                H_[14,8] = H_[8,14]
                H_[8,15] = HFPS*DP9*E8
                H_[15,8] = H_[8,15]
                H_[8,16] = HFPS*DP9*E9+DFPS*DE99
                H_[16,8] = H_[8,16]
                H_[8,17] = HFPS*DP9*ET
                H_[17,8] = H_[8,17]
                H_[9,9] = HFPS*DPT*DPT+DFPS*HPTT
                H_[9,10] = HFPS*DPT*E1
                H_[10,9] = H_[9,10]
                H_[9,11] = HFPS*DPT*E2
                H_[11,9] = H_[9,11]
                H_[9,12] = HFPS*DPT*E3
                H_[12,9] = H_[9,12]
                H_[9,13] = HFPS*DPT*E4
                H_[13,9] = H_[9,13]
                H_[9,5] = HFPS*DPT*E6
                H_[5,9] = H_[9,5]
                H_[9,14] = HFPS*DPT*E7
                H_[14,9] = H_[9,14]
                H_[9,15] = HFPS*DPT*E8
                H_[15,9] = H_[9,15]
                H_[9,16] = HFPS*DPT*E9
                H_[16,9] = H_[9,16]
                H_[9,17] = HFPS*DPT*ET+DFPS*DETT
                H_[17,9] = H_[9,17]
                H_[10,10] = HFPS*E1*E1
                H_[10,11] = HFPS*E1*E2
                H_[11,10] = H_[10,11]
                H_[10,12] = HFPS*E1*E3
                H_[12,10] = H_[10,12]
                H_[10,13] = HFPS*E1*E4
                H_[13,10] = H_[10,13]
                H_[10,5] = HFPS*E1*E6
                H_[5,10] = H_[10,5]
                H_[10,14] = HFPS*E1*E7
                H_[14,10] = H_[10,14]
                H_[10,15] = HFPS*E1*E8
                H_[15,10] = H_[10,15]
                H_[10,16] = HFPS*E1*E9
                H_[16,10] = H_[10,16]
                H_[10,17] = HFPS*E1*ET
                H_[17,10] = H_[10,17]
                H_[11,11] = HFPS*E2*E2
                H_[11,12] = HFPS*E2*E3
                H_[12,11] = H_[11,12]
                H_[11,13] = HFPS*E2*E4
                H_[13,11] = H_[11,13]
                H_[11,5] = HFPS*E2*E6
                H_[5,11] = H_[11,5]
                H_[11,14] = HFPS*E2*E7
                H_[14,11] = H_[11,14]
                H_[11,15] = HFPS*E2*E8
                H_[15,11] = H_[11,15]
                H_[11,16] = HFPS*E2*E9
                H_[16,11] = H_[11,16]
                H_[11,17] = HFPS*E2*ET
                H_[17,11] = H_[11,17]
                H_[12,12] = HFPS*E3*E3
                H_[12,13] = HFPS*E3*E4
                H_[13,12] = H_[12,13]
                H_[12,5] = HFPS*E3*E6
                H_[5,12] = H_[12,5]
                H_[12,14] = HFPS*E3*E7
                H_[14,12] = H_[12,14]
                H_[12,15] = HFPS*E3*E8
                H_[15,12] = H_[12,15]
                H_[12,16] = HFPS*E3*E9
                H_[16,12] = H_[12,16]
                H_[12,17] = HFPS*E3*ET
                H_[17,12] = H_[12,17]
                H_[13,13] = HFPS*E4*E4
                H_[13,5] = HFPS*E4*E6
                H_[5,13] = H_[13,5]
                H_[13,14] = HFPS*E4*E7
                H_[14,13] = H_[13,14]
                H_[13,15] = HFPS*E4*E8
                H_[15,13] = H_[13,15]
                H_[13,16] = HFPS*E4*E9
                H_[16,13] = H_[13,16]
                H_[13,17] = HFPS*E4*ET
                H_[17,13] = H_[13,17]
                H_[5,5] = HFPS*E6*E6
                H_[5,14] = HFPS*E6*E7
                H_[14,5] = H_[5,14]
                H_[5,15] = HFPS*E6*E8
                H_[15,5] = H_[5,15]
                H_[5,16] = HFPS*E6*E9
                H_[16,5] = H_[5,16]
                H_[5,17] = HFPS*E6*ET
                H_[17,5] = H_[5,17]
                H_[14,14] = HFPS*E7*E7
                H_[14,15] = HFPS*E7*E8
                H_[15,14] = H_[14,15]
                H_[14,16] = HFPS*E7*E9
                H_[16,14] = H_[14,16]
                H_[14,17] = HFPS*E7*ET
                H_[17,14] = H_[14,17]
                H_[15,15] = HFPS*E8*E8
                H_[15,16] = HFPS*E8*E9
                H_[16,15] = H_[15,16]
                H_[15,17] = HFPS*E8*ET
                H_[17,15] = H_[15,17]
                H_[16,16] = HFPS*E9*E9
                H_[16,17] = HFPS*E9*ET
                H_[17,16] = H_[16,17]
                H_[17,17] = HFPS*ET*ET
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

