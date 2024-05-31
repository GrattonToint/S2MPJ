from s2xlib import *
class  FCCU(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : FCCU
#    *********
# 
#    A simple data reconciliation for a fluid catalytic cracker.
# 
#                      +--------------------------+
#                      | FCCU data reconciliation |
#     +-----------+    +--------------------------+
#  1) | Flowsheet |                      Off_gas                      |------->
#     +-----------+                     |----------->                 | Propane
#                              MF_ohd   |                  |-------->F7
#                              |------>F4<-------|         |          |
#                              |        |        |         |DC3_feed  | Butane
#  Feed      Effluent          ^        |        |         |          |------->
#  ------->F1-------->F2......>|        |---------------->F5                  
#           ^                  v       DC4_feed  |         |DC4_btms  |------->
#           |                  |                 |         |          |  LCN
#           |                  |                 |<-------F6-------->F8
#           |                  | HCN              Lean_oil   C8spl_fd |  MCN
#           |                  |-------->                             |------->
#           |                  | LCO          
#           |                  |-------->              
#           |                  | HCO          
#           |                  |-------->              
#           |                  | MF_btms     
#           |                  v 
#           |<----------------F3-------->
#              Dec_recy         Decant
#     +--------------------+
#  2) | Objective function |
#     +--------------------+
#    Obj = sum 1->i [W_i*(C_flow_i - M_flow_i)**2]
#                           
#     Where: W_i       Weight on term i of objective function
#            C_flow_i  Computed flow i (a variable for this problem)
#            M_flow_i  Measrued flow i (a constant for this problem)
#     +-------------+
#  3) | Constraints |
#     +-------------+
#     These represent the linear mass balances around each
#     node, where a node (Fx) represents a single unit operation
#     in a fluid catalytics cracker.
#     +---------------+
#  4) | Initial point |
#     +---------------+
#     Feed       1.0
#     Effluent   1.0
#     MF_ohd     1.0
#     HCN        1.0
#     LCO        1.0
#     HCO        1.0
#     MF_btms    1.0
#     Decant     1.0
#     Dec_recy   1.0
#     Off_gas    1.0
#     DC4_feed   1.0
#     DC3_feed   1.0
#     DC4_btms   1.0
#     Lean_oil   1.0
#     Propane    1.0
#     Butane     1.0
#     C8spl_fd   1.0
#     LCN        1.0
#     MCN        1.0
#     Obj        7.36259000271320D+03
#     +------------------+
#  5) | Optimal solution |
#     +------------------+
#     Feed       3.11639D+01
#     Effluent   3.53528D+01
#     MF_ohd     1.94669D+01
#     HCN        2.94255D+00
#     LCO        4.94255D+00
#     HCO        3.44255D+00
#     MF_btms    4.55828D+00
#     Decant     3.69371D-01
#     Dec_recy   4.18891D+00
#     Off_gas    2.56075D+00
#     DC4_feed   2.41207D+01
#     DC3_feed   5.15601D+00
#     DC4_btms   1.89647D+01
#     Lean_oil   7.21458D+00
#     Propane    2.42801D+00
#     Butane     2.72801D+00
#     C8spl_fd   1.17501D+01
#     LCN        5.87506D+00
#     MCN        5.87506D+00
#     Obj        1.11491D+01
#     +-----------------------------------------------+
#  6) | SPEC.SPC (remove 1st * of every line to use). |
#     +-----------------------------------------------+
# BEGIN
# * maximizer-sought
# *  check-derivatives
#   ignore-derivative-bugs
# * use-scalings
# * print-scalings
#   finite-difference-gradients
# *  exact-second-derivatives-used
# * bfgs-approximate-second-derivatives-used
# * sr1-approximate-second-derivatives-used
#   bandsolver-preconditioned-cg-solver-used   5
# * diagonal-preconditioned-cg-solver-used
# * gill-murray-ponceleon-saunders-preconditioned-cg-solver-used
# * schnabel-eskow-preconditioned-cg-solver-used
# * munksgaards-preconditioned-cg-solver-used
#   exact-cauchy-point-required
# * inexact-cauchy-point-required
# * solve-bqp-accurately
# * two-norm-trust-region
# * gradient-tolerance    1.0D-5
# * constraint-tolerance  1.0D-5
#   trust-region-radius   1.0D+0
#   maximum-number-of-iterations   1000
#   print-level                    1
#   start-printing-at-iteration    0
#   stop-printing-at-iteration     1000
# END
# 
#    Source:
#    W. J. Korchinski, Profimatics, Inc,
#    325 Rolling Oaks Drive, Thousand Oaks, California, USA 91361-1200
#    Telephone: 1-805 496 6661, Fax: 1-805 373 5108
# 
#    SIF input: W. Korchinski, Spring 1993.
# 
#    classification = "SLR2-MN-19-8"
# 
# ***************************************************************
#  PROBLEM SPECIFICATION BEGINS HERE.
#  **********************************
#  **********************************           
# *************************************
#  Define objective function weights. *
# *************************************
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'FCCU'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'FCCU'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['W1'] = 0.2
        v_['W2'] = 1.0
        v_['W3'] = 1.0
        v_['W4'] = 0.33333333
        v_['W5'] = 0.33333333
        v_['W6'] = 0.33333333
        v_['W7'] = 1.0
        v_['W8'] = 1.0
        v_['W9'] = 1.0
        v_['W10'] = 1.0
        v_['W11'] = 1.0
        v_['W12'] = 1.0
        v_['W13'] = 1.0
        v_['W14'] = 1.0
        v_['W15'] = 0.33333333
        v_['W16'] = 0.33333333
        v_['W17'] = 1.0
        v_['W18'] = 0.33333333
        v_['W19'] = 0.33333333
        v_['M1'] = 31.0
        v_['M2'] = 36.0
        v_['M3'] = 20.0
        v_['M4'] = 3.0
        v_['M5'] = 5.0
        v_['M6'] = 3.5
        v_['M7'] = 4.2
        v_['M8'] = 0.9
        v_['M9'] = 3.9
        v_['M10'] = 2.2
        v_['M11'] = 22.8
        v_['M12'] = 6.8
        v_['M13'] = 19.0
        v_['M14'] = 8.5
        v_['M15'] = 2.2
        v_['M16'] = 2.5
        v_['M17'] = 10.8
        v_['M18'] = 6.5
        v_['M19'] = 6.5
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2x_ii('Feed',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Feed')
        [iv,ix_,_] = s2x_ii('Effluent',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Effluent')
        [iv,ix_,_] = s2x_ii('MFuohd',ix_)
        pb.xnames=arrset(pb.xnames,iv,'MFuohd')
        [iv,ix_,_] = s2x_ii('HCN',ix_)
        pb.xnames=arrset(pb.xnames,iv,'HCN')
        [iv,ix_,_] = s2x_ii('LCO',ix_)
        pb.xnames=arrset(pb.xnames,iv,'LCO')
        [iv,ix_,_] = s2x_ii('HCO',ix_)
        pb.xnames=arrset(pb.xnames,iv,'HCO')
        [iv,ix_,_] = s2x_ii('MFubtms',ix_)
        pb.xnames=arrset(pb.xnames,iv,'MFubtms')
        [iv,ix_,_] = s2x_ii('Decant',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Decant')
        [iv,ix_,_] = s2x_ii('Decurecy',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Decurecy')
        [iv,ix_,_] = s2x_ii('Offugas',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Offugas')
        [iv,ix_,_] = s2x_ii('DC4ufeed',ix_)
        pb.xnames=arrset(pb.xnames,iv,'DC4ufeed')
        [iv,ix_,_] = s2x_ii('DC3ufeed',ix_)
        pb.xnames=arrset(pb.xnames,iv,'DC3ufeed')
        [iv,ix_,_] = s2x_ii('DC4ubtms',ix_)
        pb.xnames=arrset(pb.xnames,iv,'DC4ubtms')
        [iv,ix_,_] = s2x_ii('Leanuoil',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Leanuoil')
        [iv,ix_,_] = s2x_ii('Propane',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Propane')
        [iv,ix_,_] = s2x_ii('Butane',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Butane')
        [iv,ix_,_] = s2x_ii('C8splufd',ix_)
        pb.xnames=arrset(pb.xnames,iv,'C8splufd')
        [iv,ix_,_] = s2x_ii('LCN',ix_)
        pb.xnames=arrset(pb.xnames,iv,'LCN')
        [iv,ix_,_] = s2x_ii('MCN',ix_)
        pb.xnames=arrset(pb.xnames,iv,'MCN')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2x_ii('F1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'F1')
        iv = ix_['Feed']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['Decurecy']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['Effluent']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('F2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'F2')
        iv = ix_['Effluent']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['MFuohd']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['HCN']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['LCO']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['HCO']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['MFubtms']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('F3',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'F3')
        iv = ix_['MFubtms']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['Decant']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['Decurecy']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('F4',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'F4')
        iv = ix_['MFuohd']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['Leanuoil']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['Offugas']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['DC4ufeed']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('F5',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'F5')
        iv = ix_['DC4ufeed']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['DC3ufeed']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['DC4ubtms']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('F6',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'F6')
        iv = ix_['DC4ubtms']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['Leanuoil']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['C8splufd']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('F7',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'F7')
        iv = ix_['DC3ufeed']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['Propane']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['Butane']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('F8',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'F8')
        iv = ix_['C8splufd']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['LCN']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['MCN']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('Obj1',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['Feed']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['W1']))
        [ig,ig_,_] = s2x_ii('Obj2',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['Effluent']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['W2']))
        [ig,ig_,_] = s2x_ii('Obj3',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['MFuohd']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['W3']))
        [ig,ig_,_] = s2x_ii('Obj4',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['HCN']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['W4']))
        [ig,ig_,_] = s2x_ii('Obj5',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['LCO']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['W5']))
        [ig,ig_,_] = s2x_ii('Obj6',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['HCO']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['W6']))
        [ig,ig_,_] = s2x_ii('Obj7',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['MFubtms']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['W7']))
        [ig,ig_,_] = s2x_ii('Obj8',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['Decant']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['W8']))
        [ig,ig_,_] = s2x_ii('Obj9',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['Decurecy']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['W9']))
        [ig,ig_,_] = s2x_ii('Obj10',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['Offugas']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['W10']))
        [ig,ig_,_] = s2x_ii('Obj11',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['DC4ufeed']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['W11']))
        [ig,ig_,_] = s2x_ii('Obj12',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['DC3ufeed']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['W12']))
        [ig,ig_,_] = s2x_ii('Obj13',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['DC4ubtms']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['W13']))
        [ig,ig_,_] = s2x_ii('Obj14',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['Leanuoil']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['W14']))
        [ig,ig_,_] = s2x_ii('Obj15',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['Propane']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['W15']))
        [ig,ig_,_] = s2x_ii('Obj16',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['Butane']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['W16']))
        [ig,ig_,_] = s2x_ii('Obj17',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['C8splufd']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['W17']))
        [ig,ig_,_] = s2x_ii('Obj18',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['LCN']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['W18']))
        [ig,ig_,_] = s2x_ii('Obj19',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['MCN']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['W19']))
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
        pbm.gconst = arrset(pbm.gconst,ig_['Obj1'],float(v_['M1']))
        pbm.gconst = arrset(pbm.gconst,ig_['Obj2'],float(v_['M2']))
        pbm.gconst = arrset(pbm.gconst,ig_['Obj3'],float(v_['M3']))
        pbm.gconst = arrset(pbm.gconst,ig_['Obj4'],float(v_['M4']))
        pbm.gconst = arrset(pbm.gconst,ig_['Obj5'],float(v_['M5']))
        pbm.gconst = arrset(pbm.gconst,ig_['Obj6'],float(v_['M6']))
        pbm.gconst = arrset(pbm.gconst,ig_['Obj7'],float(v_['M7']))
        pbm.gconst = arrset(pbm.gconst,ig_['Obj8'],float(v_['M8']))
        pbm.gconst = arrset(pbm.gconst,ig_['Obj9'],float(v_['M9']))
        pbm.gconst = arrset(pbm.gconst,ig_['Obj10'],float(v_['M10']))
        pbm.gconst = arrset(pbm.gconst,ig_['Obj11'],float(v_['M11']))
        pbm.gconst = arrset(pbm.gconst,ig_['Obj12'],float(v_['M12']))
        pbm.gconst = arrset(pbm.gconst,ig_['Obj13'],float(v_['M13']))
        pbm.gconst = arrset(pbm.gconst,ig_['Obj14'],float(v_['M14']))
        pbm.gconst = arrset(pbm.gconst,ig_['Obj15'],float(v_['M15']))
        pbm.gconst = arrset(pbm.gconst,ig_['Obj16'],float(v_['M16']))
        pbm.gconst = arrset(pbm.gconst,ig_['Obj17'],float(v_['M17']))
        pbm.gconst = arrset(pbm.gconst,ig_['Obj18'],float(v_['M18']))
        pbm.gconst = arrset(pbm.gconst,ig_['Obj19'],float(v_['M19']))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        pb.x0[ix_['Feed']] = float(1.0)
        pb.x0[ix_['Effluent']] = float(1.0)
        pb.x0[ix_['MFuohd']] = float(1.0)
        pb.x0[ix_['HCN']] = float(1.0)
        pb.x0[ix_['LCO']] = float(1.0)
        pb.x0[ix_['HCO']] = float(1.0)
        pb.x0[ix_['MFubtms']] = float(1.0)
        pb.x0[ix_['Decant']] = float(1.0)
        pb.x0[ix_['Decurecy']] = float(1.0)
        pb.x0[ix_['Offugas']] = float(1.0)
        pb.x0[ix_['DC4ufeed']] = float(1.0)
        pb.x0[ix_['DC3ufeed']] = float(1.0)
        pb.x0[ix_['DC4ubtms']] = float(1.0)
        pb.x0[ix_['Leanuoil']] = float(1.0)
        pb.x0[ix_['Propane']] = float(1.0)
        pb.x0[ix_['Butane']] = float(1.0)
        pb.x0[ix_['C8splufd']] = float(1.0)
        pb.x0[ix_['LCN']] = float(1.0)
        pb.x0[ix_['MCN']] = float(1.0)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2x_ii('gSQUARE',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        ig = ig_['Obj1']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        ig = ig_['Obj2']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        ig = ig_['Obj3']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        ig = ig_['Obj4']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        ig = ig_['Obj5']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        ig = ig_['Obj6']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        ig = ig_['Obj7']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        ig = ig_['Obj8']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        ig = ig_['Obj9']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        ig = ig_['Obj10']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        ig = ig_['Obj11']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        ig = ig_['Obj12']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        ig = ig_['Obj13']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        ig = ig_['Obj14']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        ig = ig_['Obj15']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        ig = ig_['Obj16']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        ig = ig_['Obj17']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        ig = ig_['Obj18']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        ig = ig_['Obj19']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
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
        pb.pbclass = "SLR2-MN-19-8"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gSQUARE(pbm,nargout,*args):

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
