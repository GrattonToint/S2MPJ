from s2mpjlib import *
class  LINSPANH(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : LINSPANH
#    *********
# 
#    A linear network problem based on the spanish hydro-electric
#    reservoir management problem SPANHYD
# 
#    Source:
#    A partial specification of problem SPANHYD.
# 
#    SIF input: Ph. Toint, Sept 1990.
# 
#    classification = "C-CLNR2-MN-97-33"
# 
#    Number of arcs = 97
#    Number of nodes
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'LINSPANH'

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
        irA  = np.append(irA,[ig_['OBJ']])
        valA = np.append(valA,float(-1.0))
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
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN                77.0
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons   = np.arange(len(self.congrps))
        self.pbclass   = "C-CLNR2-MN-97-33"
        self.objderlvl = 2
        self.conderlvl = [2]


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

