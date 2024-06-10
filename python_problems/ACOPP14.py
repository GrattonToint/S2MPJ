from s2mpjlib import *
class  ACOPP14(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : ACOPP14
#    *********
# 
#    An AC Optimal Power Flow (OPF) problem for the IEEE 14 Bus
#    Power Systems Test Case from the archive:
#      http://www.ee.washington.edu/research/pstca/
# 
#    Polar formulation due to 
#     Anya Castillo, Johns Hopkins University, anya.castillo@jhu.edu
# 
#    variables: 
#      A - voltage amplitude
#      M (= |V|) - voltage modulus
#      P - real power component
#      Q - imaginary (reactive) power component
# 
#    constants  
#      PD, QD - constant loads (= withdrawl of power from the network)
#      Pmin, Pmax - real power limits
#      Qmin, Qmax - imaginary realtive) power limits
#      Mmin, Mmax - voltage modulus limits
#      Amin, Amax - voltage phase aplitude limits
#      Smax - thermal limits
# 
#    objective function:
#    ------------------
#      sum_{i in nodes} f_i(P_i) = a_i P_i^2 + b_i P_i
# 
#    real power flow constraints:
#    ---------------------------
#      R_i * ( G * R - B * I ) + I_i * ( G * I + B * R ) 
#        - P_i + PD_i = 0 for all nodes i
#      M_i * sum_{j in nodes} M_j * ( G_ij cos A_ij  + B_ij sin A_ij ) 
#        - P_i + PD_i = 0 for all nodes i
#      where A_ij = A_i - A_j
# 
#    reactive power flow constraints:
#    -------------------------------
#      I_i * ( G * R - B * I ) - R_i * ( G * I + B * R ) 
#        - Q_i + QD_i = 0 for all nodes i
#      M_i * sum_{j in nodes} M_j * ( G_ij sin A_ij  - B_ij cos A_ij ) 
#        - Q_i + QD_i = 0 for all nodes i
# 
#    line thermal limit constraints:
#    ------------------------------
#      f_i(A,M) <= Smax_i and  t_i(A,M) <= Smax_i for all lines i
# 
#      here if we write v = M * ( cos A + i sin A ) = v^R + i V^I,
#        f_i = ( v_j(i) . ( Yf * v )_i ) * conj( v_j(i) . ( Yf * v )_i )
#        t_i = ( v_j(i) . ( Yt * v )_i ) * conj( v_j(i) . ( Yt * v )_i )
#      where the nodes j(i) and 
#        Yf = Yf^R + i Yf^I and Yt = Yt^R + i Yt^I are given
#      This leads to
#         f_i = ( R_j . Yf^R R - R_j . Yf^I I + I_j . Yf^R I + I_j . Yf^I R )^2 +
#               ( I_j . Yf^R R - I_j . Yf^I I - R_j . Yf^R I - R_j . Yf^I R )^2
#       and
#         t_i = ( R_j . Yt^R R - R_j . Yt^I I + I_j . Yt^R I + I_j . Yt^I R )^2 +
#               ( I_j . Yt^R R - I_j . Yt^I I - R_j . Yt^R I - R_j . Yt^I R )^2
# 
#      if Yt^R R = sum_k Yt_k^R R_k (etc), we have
# 
#         f_i = ( sum_k [ R_j . Yf^R_k R_k - R_j . Yf^I_k I_k + 
#                         I_j . Yf^R_k I_k + I_j . Yf^I_k R_k ] )^2 +
#               ( sum_k [ I_j . Yf^R_k R_k - I_j . Yf^I_k I_k - 
#                         R_j . Yf^R_k I_k - R_j . Yf^I_k R_k ] )^2
#             =  sum_k [          ( Yf^R_k^2 + Yf^I_k^2 ) R_j^2 R_k^2
#                               + ( Yf^R_k^2 + Yf^I_k^2 ) I_j^2 I_k^2
#                               + ( Yf^R_k^2 + Yf^I_k^2 ) R_j^2 I_k^2
#                               + ( Yf^R_k^2 + Yf^I_k^2 ) I_j^2 R_k^2 ] +
#               sum_k sum_l>k [   2 ( Yf^R_k Yf^R_l + Yf^I_k Yf^I_l ) 
#                                     R_j^2 R_k R_l 
#                               + 2 ( Yf^R_k Yf^R_l + Yf^I_k Yf^I_l ) 
#                                     R_j^2 I_k I_l 
#                               + 2 ( Yf^R_k Yf^R_l + Yf^I_k Yf^I_l ) 
#                                     I_j^2 R_k R_l 
#                               + 2 ( Yf^R_k Yf^R_l + Yf^I_k Yf^I_l ) 
#                                     I_j^2 I_k I_l 
#                               + 2 ( Yf^R_k Yf^I_l - Yf^I_k Yf^R_l )
#                                     R_j^2 I_k R_l 
#                               - 2 ( Yf^R_k Yf^I_l - Yf^I_k Yf^R_l )
#                                     R_j^2 R_k I_l 
#                               + 2 ( Yf^R_k Yf^I_l - Yf^I_k Yf^R_l )
#                                     I_j^2 I_k R_l 
#                               - 2 ( Yf^R_k Yf^I_l - Yf^I_k Yf^R_l )
#                                     I_j^2 R_k I_l                   ]
#             = sum_k                ( Yf^R_k^2 + Yf^I_k^2 ) M_j^2 M_k^2 +
#               sum_k sum_l>k [    2 ( Yf^R_k Yf^R_l + Yf^I_k Yf^I_l ) 
#                                     M_j^2 M_k M_l cos ( A_k - A_l )
#                               +  2 ( Yf^R_k Yf^I_l - Yf^I_k Yf^R_l ) 
#                                     M_j^2 M_k M_l sin ( A_k - A_l ) ]
# 
#      and similarly for t_i
# 
#  [ ** NOT USED **
#    maximum phase-amplitude difference constraints:
#      Amin_ij <= A_i - A_j <= Amax_ij  for all interconnets i and j ] 
# 
#    node voltage modulus limits:
#    ----------------------------
#      Mmin_i <= M_i <= Mmax_i  for all nodes i
# 
#    generator real power limits:
#    ---------------------------
#      Pmin_i <= P_i <= Pmax_i  for all nodes i
# 
#    generator reactive power limits:
#    -------------------------------
#      Qmin_i <= Q_i <= Qmax_i  for all nodes i
# 
#    SIF input: Nick Gould, August 2011
# 
#    classification = "QOR2-AY-38-68"
# 
#    number of nodes
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'ACOPP14'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['NODES'] = 14
        v_['LIMITS'] = 5
        v_['LINES'] = 20
        v_['1'] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['NODES'])+1):
            [iv,ix_,_] = s2mpj_ii('A'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'A'+str(I))
            [iv,ix_,_] = s2mpj_ii('M'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'M'+str(I))
        for I in range(int(v_['1']),int(v_['LIMITS'])+1):
            [iv,ix_,_] = s2mpj_ii('P'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'P'+str(I))
            [iv,ix_,_] = s2mpj_ii('Q'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'Q'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['P1']
        pbm.A[ig,iv] = float(2000.0)+pbm.A[ig,iv]
        iv = ix_['P2']
        pbm.A[ig,iv] = float(2000.0)+pbm.A[ig,iv]
        iv = ix_['P3']
        pbm.A[ig,iv] = float(4000.0)+pbm.A[ig,iv]
        iv = ix_['P4']
        pbm.A[ig,iv] = float(4000.0)+pbm.A[ig,iv]
        iv = ix_['P5']
        pbm.A[ig,iv] = float(4000.0)+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['NODES'])+1):
            [ig,ig_,_] = s2mpj_ii('RP'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'RP'+str(I))
        for I in range(int(v_['1']),int(v_['NODES'])+1):
            [ig,ig_,_] = s2mpj_ii('IP'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'IP'+str(I))
        [ig,ig_,_] = s2mpj_ii('RP1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'RP1')
        iv = ix_['P1']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('IP1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'IP1')
        iv = ix_['Q1']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('RP2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'RP2')
        iv = ix_['P2']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('IP2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'IP2')
        iv = ix_['Q2']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('RP3',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'RP3')
        iv = ix_['P3']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('IP3',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'IP3')
        iv = ix_['Q3']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('RP6',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'RP6')
        iv = ix_['P4']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('IP6',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'IP6')
        iv = ix_['Q4']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('RP8',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'RP8')
        iv = ix_['P5']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('IP8',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'IP8')
        iv = ix_['Q5']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['LINES'])+1):
            [ig,ig_,_] = s2mpj_ii('FN'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'FN'+str(I))
        for I in range(int(v_['1']),int(v_['LINES'])+1):
            [ig,ig_,_] = s2mpj_ii('TN'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'TN'+str(I))
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
        pbm.gconst = arrset(pbm.gconst,ig_['RP2'],float(-0.217))
        pbm.gconst = arrset(pbm.gconst,ig_['RP3'],float(-0.942))
        pbm.gconst = arrset(pbm.gconst,ig_['RP4'],float(-0.478))
        pbm.gconst = arrset(pbm.gconst,ig_['RP5'],float(-0.076))
        pbm.gconst = arrset(pbm.gconst,ig_['RP6'],float(-0.112))
        pbm.gconst = arrset(pbm.gconst,ig_['RP9'],float(-0.295))
        pbm.gconst = arrset(pbm.gconst,ig_['RP10'],float(-0.09))
        pbm.gconst = arrset(pbm.gconst,ig_['RP11'],float(-0.035))
        pbm.gconst = arrset(pbm.gconst,ig_['RP12'],float(-0.061))
        pbm.gconst = arrset(pbm.gconst,ig_['RP13'],float(-0.135))
        pbm.gconst = arrset(pbm.gconst,ig_['RP14'],float(-0.149))
        pbm.gconst = arrset(pbm.gconst,ig_['IP2'],float(-0.127))
        pbm.gconst = arrset(pbm.gconst,ig_['IP3'],float(-0.19))
        pbm.gconst = arrset(pbm.gconst,ig_['IP4'],float(0.039))
        pbm.gconst = arrset(pbm.gconst,ig_['IP5'],float(-0.016))
        pbm.gconst = arrset(pbm.gconst,ig_['IP6'],float(-0.075))
        pbm.gconst = arrset(pbm.gconst,ig_['IP9'],float(-0.166))
        pbm.gconst = arrset(pbm.gconst,ig_['IP10'],float(-0.058))
        pbm.gconst = arrset(pbm.gconst,ig_['IP11'],float(-0.018))
        pbm.gconst = arrset(pbm.gconst,ig_['IP12'],float(-0.016))
        pbm.gconst = arrset(pbm.gconst,ig_['IP13'],float(-0.058))
        pbm.gconst = arrset(pbm.gconst,ig_['IP14'],float(-0.05))
        pbm.gconst = arrset(pbm.gconst,ig_['FN1'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['FN2'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['FN3'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['FN4'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['FN5'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['FN6'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['FN7'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['FN8'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['FN9'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['FN10'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['FN11'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['FN12'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['FN13'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['FN14'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['FN15'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['FN16'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['FN17'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['FN18'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['FN19'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['FN20'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['TN1'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['TN2'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['TN3'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['TN4'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['TN5'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['TN6'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['TN7'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['TN8'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['TN9'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['TN10'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['TN11'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['TN12'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['TN13'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['TN14'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['TN15'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['TN16'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['TN17'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['TN18'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['TN19'],float(9801.0))
        pbm.gconst = arrset(pbm.gconst,ig_['TN20'],float(9801.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-1.0e+30)
        pb.xupper = np.full((pb.n,1),1.0e+30)
        pb.xlower[ix_['M1']] = 0.94
        pb.xupper[ix_['M1']] = 1.06
        pb.xlower[ix_['M2']] = 0.94
        pb.xupper[ix_['M2']] = 1.06
        pb.xlower[ix_['M3']] = 0.94
        pb.xupper[ix_['M3']] = 1.06
        pb.xlower[ix_['M4']] = 0.94
        pb.xupper[ix_['M4']] = 1.06
        pb.xlower[ix_['M5']] = 0.94
        pb.xupper[ix_['M5']] = 1.06
        pb.xlower[ix_['M6']] = 0.94
        pb.xupper[ix_['M6']] = 1.06
        pb.xlower[ix_['M7']] = 0.94
        pb.xupper[ix_['M7']] = 1.06
        pb.xlower[ix_['M8']] = 0.94
        pb.xupper[ix_['M8']] = 1.06
        pb.xlower[ix_['M9']] = 0.94
        pb.xupper[ix_['M9']] = 1.06
        pb.xlower[ix_['M10']] = 0.94
        pb.xupper[ix_['M10']] = 1.06
        pb.xlower[ix_['M11']] = 0.94
        pb.xupper[ix_['M11']] = 1.06
        pb.xlower[ix_['M12']] = 0.94
        pb.xupper[ix_['M12']] = 1.06
        pb.xlower[ix_['M13']] = 0.94
        pb.xupper[ix_['M13']] = 1.06
        pb.xlower[ix_['M14']] = 0.94
        pb.xupper[ix_['M14']] = 1.06
        pb.xlower[ix_['P1']] = 0.0
        pb.xupper[ix_['P1']] = 3.324
        pb.xlower[ix_['P2']] = 0.0
        pb.xupper[ix_['P2']] = 1.4
        pb.xlower[ix_['P3']] = 0.0
        pb.xupper[ix_['P3']] = 1.0
        pb.xlower[ix_['P4']] = 0.0
        pb.xupper[ix_['P4']] = 1.0
        pb.xlower[ix_['P5']] = 0.0
        pb.xupper[ix_['P5']] = 1.0
        pb.xlower[ix_['Q1']] = 0.0
        pb.xupper[ix_['Q1']] = 0.1
        pb.xlower[ix_['Q2']] = -0.4
        pb.xupper[ix_['Q2']] = 0.5
        pb.xlower[ix_['Q3']] = 0.0
        pb.xupper[ix_['Q3']] = 0.4
        pb.xlower[ix_['Q4']] = -0.06
        pb.xupper[ix_['Q4']] = 0.24
        pb.xlower[ix_['Q5']] = -0.06
        pb.xupper[ix_['Q5']] = 0.24
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(0.0))
        pb.x0[ix_['M1']] = float(1.06)
        pb.x0[ix_['M2']] = float(1.045)
        pb.x0[ix_['M3']] = float(1.01)
        pb.x0[ix_['M4']] = float(1.019)
        pb.x0[ix_['M5']] = float(1.02)
        pb.x0[ix_['M6']] = float(1.07)
        pb.x0[ix_['M7']] = float(1.062)
        pb.x0[ix_['M8']] = float(1.09)
        pb.x0[ix_['M9']] = float(1.056)
        pb.x0[ix_['M10']] = float(1.051)
        pb.x0[ix_['M11']] = float(1.057)
        pb.x0[ix_['M12']] = float(1.055)
        pb.x0[ix_['M13']] = float(1.05)
        pb.x0[ix_['M14']] = float(1.036)
        pb.x0[ix_['A2']] = float(-0.086917396)
        pb.x0[ix_['A3']] = float(-0.222005880)
        pb.x0[ix_['A4']] = float(-0.180292511)
        pb.x0[ix_['A5']] = float(-0.153239908)
        pb.x0[ix_['A6']] = float(-0.248185819)
        pb.x0[ix_['A7']] = float(-0.233350520)
        pb.x0[ix_['A8']] = float(-0.233175988)
        pb.x0[ix_['A9']] = float(-0.260752190)
        pb.x0[ix_['A10']] = float(-0.263544717)
        pb.x0[ix_['A11']] = float(-0.258134196)
        pb.x0[ix_['A12']] = float(-0.263021118)
        pb.x0[ix_['A13']] = float(-0.264591914)
        pb.x0[ix_['A14']] = float(-0.279950812)
        pb.x0[ix_['P1']] = float(2.324)
        pb.x0[ix_['P2']] = float(0.4)
        pb.x0[ix_['Q1']] = float(-0.169)
        pb.x0[ix_['Q2']] = float(0.424)
        pb.x0[ix_['Q3']] = float(0.234)
        pb.x0[ix_['Q4']] = float(0.122)
        pb.x0[ix_['Q5']] = float(0.174)
        #%%%%%%%%%%%%%%%%%%%% QUADRATIC %%%%%%%%%%%%%%%%%%%
        pbm.H = lil_matrix((pb.n, pb.n))
        ix1 = ix_['P1']
        ix2 = ix_['P1']
        pbm.H[ix1,ix2] = float(860.586)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_['P2']
        ix2 = ix_['P2']
        pbm.H[ix1,ix2] = float(5000.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_['P3']
        ix2 = ix_['P3']
        pbm.H[ix1,ix2] = float(200.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_['P4']
        ix2 = ix_['P4']
        pbm.H[ix1,ix2] = float(200.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        ix1 = ix_['P5']
        ix2 = ix_['P5']
        pbm.H[ix1,ix2] = float(200.0)+pbm.H[ix1,ix2]
        pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eP2', iet_)
        elftv = loaset(elftv,it,0,'V1')
        [it,iet_,_] = s2mpj_ii( 'eP4', iet_)
        elftv = loaset(elftv,it,0,'V1')
        [it,iet_,_] = s2mpj_ii( 'eP11', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        [it,iet_,_] = s2mpj_ii( 'eP31', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        [it,iet_,_] = s2mpj_ii( 'eP22', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        [it,iet_,_] = s2mpj_ii( 'eP211', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        [it,iet_,_] = s2mpj_ii( 'eSIN11', iet_)
        elftv = loaset(elftv,it,0,'U1')
        elftv = loaset(elftv,it,1,'U2')
        elftv = loaset(elftv,it,2,'A1')
        elftv = loaset(elftv,it,3,'A2')
        [it,iet_,_] = s2mpj_ii( 'eCOS11', iet_)
        elftv = loaset(elftv,it,0,'U1')
        elftv = loaset(elftv,it,1,'U2')
        elftv = loaset(elftv,it,2,'A1')
        elftv = loaset(elftv,it,3,'A2')
        [it,iet_,_] = s2mpj_ii( 'eSIN211', iet_)
        elftv = loaset(elftv,it,0,'U1')
        elftv = loaset(elftv,it,1,'U2')
        elftv = loaset(elftv,it,2,'U3')
        elftv = loaset(elftv,it,3,'A1')
        elftv = loaset(elftv,it,4,'A2')
        [it,iet_,_] = s2mpj_ii( 'eCOS211', iet_)
        elftv = loaset(elftv,it,0,'U1')
        elftv = loaset(elftv,it,1,'U2')
        elftv = loaset(elftv,it,2,'U3')
        elftv = loaset(elftv,it,3,'A1')
        elftv = loaset(elftv,it,4,'A2')
        [it,iet_,_] = s2mpj_ii( 'eSIN31', iet_)
        elftv = loaset(elftv,it,0,'U1')
        elftv = loaset(elftv,it,1,'U2')
        elftv = loaset(elftv,it,2,'A1')
        elftv = loaset(elftv,it,3,'A2')
        [it,iet_,_] = s2mpj_ii( 'eCOS31', iet_)
        elftv = loaset(elftv,it,0,'U1')
        elftv = loaset(elftv,it,1,'U2')
        elftv = loaset(elftv,it,2,'A1')
        elftv = loaset(elftv,it,3,'A2')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        ename = 'F1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP2')
        ielftype = arrset(ielftype, ie, iet_["eP2"])
        vname = 'M1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP2')
        ielftype = arrset(ielftype, ie, iet_["eP2"])
        vname = 'M1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F5'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F6'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F7'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F8'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F9'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F10'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F11'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F12'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F13'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F14'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F15'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP2')
        ielftype = arrset(ielftype, ie, iet_["eP2"])
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F16'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP2')
        ielftype = arrset(ielftype, ie, iet_["eP2"])
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F17'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F18'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F19'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F20'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F21'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F22'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F23'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F24'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F25'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F26'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F27'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F28'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F29'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F30'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F31'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F32'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F33'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP2')
        ielftype = arrset(ielftype, ie, iet_["eP2"])
        vname = 'M3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F34'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP2')
        ielftype = arrset(ielftype, ie, iet_["eP2"])
        vname = 'M3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F35'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F36'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F37'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F38'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F39'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F40'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F41'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F42'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F43'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F44'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F45'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F46'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F47'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP2')
        ielftype = arrset(ielftype, ie, iet_["eP2"])
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F48'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP2')
        ielftype = arrset(ielftype, ie, iet_["eP2"])
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F49'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F50'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F51'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F52'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F53'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F54'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F55'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F56'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F57'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F58'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F59'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F60'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F61'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F62'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F63'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F64'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F65'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F66'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F67'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F68'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F69'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP2')
        ielftype = arrset(ielftype, ie, iet_["eP2"])
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F70'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP2')
        ielftype = arrset(ielftype, ie, iet_["eP2"])
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F71'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F72'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F73'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F74'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F75'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP2')
        ielftype = arrset(ielftype, ie, iet_["eP2"])
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F76'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP2')
        ielftype = arrset(ielftype, ie, iet_["eP2"])
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F77'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F78'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F79'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F80'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F81'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F82'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F83'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F84'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F85'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F86'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F87'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F88'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F89'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F90'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F91'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP2')
        ielftype = arrset(ielftype, ie, iet_["eP2"])
        vname = 'M7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F92'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M8'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A8'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F93'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M8'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A8'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F94'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F95'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F96'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M8'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A8'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F97'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M8'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A8'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F98'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP2')
        ielftype = arrset(ielftype, ie, iet_["eP2"])
        vname = 'M8'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F99'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F100'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F101'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F102'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F103'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP2')
        ielftype = arrset(ielftype, ie, iet_["eP2"])
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F104'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP2')
        ielftype = arrset(ielftype, ie, iet_["eP2"])
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F105'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F106'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F107'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F108'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F109'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F110'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F111'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F112'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F113'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F114'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F115'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F116'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F117'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP2')
        ielftype = arrset(ielftype, ie, iet_["eP2"])
        vname = 'M10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F118'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP2')
        ielftype = arrset(ielftype, ie, iet_["eP2"])
        vname = 'M10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F119'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F120'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F121'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F122'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F123'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F124'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F125'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F126'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F127'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F128'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F129'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F130'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F131'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP2')
        ielftype = arrset(ielftype, ie, iet_["eP2"])
        vname = 'M11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F132'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP2')
        ielftype = arrset(ielftype, ie, iet_["eP2"])
        vname = 'M11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F133'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F134'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F135'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F136'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F137'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP2')
        ielftype = arrset(ielftype, ie, iet_["eP2"])
        vname = 'M12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F138'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP2')
        ielftype = arrset(ielftype, ie, iet_["eP2"])
        vname = 'M12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F139'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F140'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F141'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F142'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F143'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F144'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F145'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F146'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F147'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F148'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F149'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F150'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F151'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP2')
        ielftype = arrset(ielftype, ie, iet_["eP2"])
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F152'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP2')
        ielftype = arrset(ielftype, ie, iet_["eP2"])
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F153'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F154'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F155'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F156'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F157'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F158'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F159'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F160'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F161'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F162'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F163'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype, ie, iet_["eSIN11"])
        vname = 'M14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F164'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype, ie, iet_["eCOS11"])
        vname = 'M14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F165'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP2')
        ielftype = arrset(ielftype, ie, iet_["eP2"])
        vname = 'M14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'F166'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP2')
        ielftype = arrset(ielftype, ie, iet_["eP2"])
        vname = 'M14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN31')
        ielftype = arrset(ielftype, ie, iet_["eSIN31"])
        vname = 'M1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E5'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E6'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E7'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN31')
        ielftype = arrset(ielftype, ie, iet_["eSIN31"])
        vname = 'M1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E8'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E9'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E10'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E11'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN31')
        ielftype = arrset(ielftype, ie, iet_["eSIN31"])
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E12'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E13'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E14'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E15'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN31')
        ielftype = arrset(ielftype, ie, iet_["eSIN31"])
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E16'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E17'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E18'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E19'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN31')
        ielftype = arrset(ielftype, ie, iet_["eSIN31"])
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E20'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E21'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E22'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E23'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN31')
        ielftype = arrset(ielftype, ie, iet_["eSIN31"])
        vname = 'M3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E24'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E25'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E26'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E27'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E28'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E29'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E30'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E31'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E32'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E33'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E34'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E35'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E36'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E37'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E38'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E39'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E40'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E41'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E42'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E43'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E44'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E45'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E46'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E47'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M8'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A8'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E48'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M8'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E49'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E50'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E51'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E52'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E53'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E54'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E55'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E56'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E57'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E58'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E59'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E60'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E61'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E62'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E63'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E64'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E65'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E66'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E67'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E68'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E69'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN31')
        ielftype = arrset(ielftype, ie, iet_["eSIN31"])
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E70'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E71'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E72'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E73'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN31')
        ielftype = arrset(ielftype, ie, iet_["eSIN31"])
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E74'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E75'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E76'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E77'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN31')
        ielftype = arrset(ielftype, ie, iet_["eSIN31"])
        vname = 'M3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E78'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E79'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E80'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E81'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN31')
        ielftype = arrset(ielftype, ie, iet_["eSIN31"])
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E82'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E83'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E84'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E85'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN31')
        ielftype = arrset(ielftype, ie, iet_["eSIN31"])
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E86'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E87'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E88'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E89'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN31')
        ielftype = arrset(ielftype, ie, iet_["eSIN31"])
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E90'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E91'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E92'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E93'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E94'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E95'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E96'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E97'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E98'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E99'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E100'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E101'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E102'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E103'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E104'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E105'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E106'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E107'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E108'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E109'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E110'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E111'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E112'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M8'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E113'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M8'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A8'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E114'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M8'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E115'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E116'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E117'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E118'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E119'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E120'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E121'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E122'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E123'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E124'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E125'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E126'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E127'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E128'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E129'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E130'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP22')
        ielftype = arrset(ielftype, ie, iet_["eP22"])
        vname = 'M14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E131'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype, ie, iet_["eCOS31"])
        vname = 'M14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'M13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A13'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E132'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eP4')
        ielftype = arrset(ielftype, ie, iet_["eP4"])
        vname = 'M14'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-1.0e+30,1.0e+30,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        ig = ig_['RP1']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F1'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(6.0250290558))
        ig = ig_['IP1']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F2'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(19.447070206))
        ig = ig_['RP1']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F3'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-4.999131600))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F4'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(15.263086523))
        ig = ig_['IP1']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F5'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-4.999131600))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F6'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-15.26308652))
        ig = ig_['RP1']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F7'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.025897455))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F8'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(4.2349836823))
        ig = ig_['IP1']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F9'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.025897455))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F10'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-4.234983682))
        ig = ig_['RP2']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F11'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-4.999131600))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F12'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(15.263086523))
        ig = ig_['IP2']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F13'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-4.999131600))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F14'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-15.26308652))
        ig = ig_['RP2']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F15'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(9.5213236108))
        ig = ig_['IP2']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F16'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(30.272115399))
        ig = ig_['RP2']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F17'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.135019192))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F18'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(4.7818631518))
        ig = ig_['IP2']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F19'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.135019192))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F20'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-4.781863151))
        ig = ig_['RP2']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F21'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.686033150))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F22'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(5.1158383259))
        ig = ig_['IP2']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F23'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.686033150))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F24'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-5.115838325))
        ig = ig_['RP2']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F25'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.701139667))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F26'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(5.193927398))
        ig = ig_['IP2']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F27'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.701139667))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F28'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-5.193927398))
        ig = ig_['RP3']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F29'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.135019192))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F30'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(4.7818631518))
        ig = ig_['IP3']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F31'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.135019192))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F32'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-4.781863151))
        ig = ig_['RP3']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F33'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(3.1209949022))
        ig = ig_['IP3']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F34'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(9.8223801294))
        ig = ig_['RP3']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F35'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.985975709))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F36'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(5.0688169776))
        ig = ig_['IP3']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F37'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.985975709))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F38'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-5.068816977))
        ig = ig_['RP4']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F39'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.686033150))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F40'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(5.1158383259))
        ig = ig_['IP4']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F41'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.686033150))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F42'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-5.115838325))
        ig = ig_['RP4']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F43'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.985975709))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F44'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(5.0688169776))
        ig = ig_['IP4']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F45'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.985975709))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F46'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-5.068816977))
        ig = ig_['RP4']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F47'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(10.512989522))
        ig = ig_['IP4']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F48'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(38.654171208))
        ig = ig_['RP4']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F49'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-6.840980661))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F50'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(21.578553982))
        ig = ig_['IP4']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F51'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-6.840980661))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F52'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-21.57855398))
        ig = ig_['RP4']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F53'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(4.8895126603))
        ig = ig_['IP4']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F54'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-4.889512660))
        ig = ig_['RP4']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F55'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.8554995578))
        ig = ig_['IP4']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F56'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.855499557))
        ig = ig_['RP5']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F57'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.025897455))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F58'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(4.2349836823))
        ig = ig_['IP5']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F59'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.025897455))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F60'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-4.234983682))
        ig = ig_['RP5']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F61'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.701139667))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F62'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(5.193927398))
        ig = ig_['IP5']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F63'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.701139667))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F64'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-5.193927398))
        ig = ig_['RP5']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F65'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-6.840980661))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F66'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(21.578553982))
        ig = ig_['IP5']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F67'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-6.840980661))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F68'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-21.57855398))
        ig = ig_['RP5']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F69'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(9.5680177836))
        ig = ig_['IP5']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F70'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(35.533639456))
        ig = ig_['RP5']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F71'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(4.2574453353))
        ig = ig_['IP5']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F72'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-4.257445335))
        ig = ig_['RP6']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F73'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(4.2574453353))
        ig = ig_['IP6']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F74'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-4.257445335))
        ig = ig_['RP6']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F75'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(6.5799234075))
        ig = ig_['IP6']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F76'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(17.34073281))
        ig = ig_['RP6']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F77'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.955028563))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F78'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(4.0940743442))
        ig = ig_['IP6']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F79'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.955028563))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F80'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-4.094074344))
        ig = ig_['RP6']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F81'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.525967440))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F82'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(3.175963965))
        ig = ig_['IP6']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F83'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.525967440))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F84'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-3.175963965))
        ig = ig_['RP6']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F85'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-3.098927403))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F86'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(6.1027554482))
        ig = ig_['IP6']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F87'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-3.098927403))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F88'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-6.102755448))
        ig = ig_['RP7']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F89'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(4.8895126603))
        ig = ig_['IP7']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F90'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-4.889512660))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F91'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(19.549005948))
        ig = ig_['RP7']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F92'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(5.6769798467))
        ig = ig_['IP7']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F93'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-5.676979846))
        ig = ig_['RP7']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F94'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(9.0900827198))
        ig = ig_['IP7']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F95'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-9.090082719))
        ig = ig_['RP8']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F96'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(5.6769798467))
        ig = ig_['IP8']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F97'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-5.676979846))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F98'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(5.6769798467))
        ig = ig_['RP9']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F99'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.8554995578))
        ig = ig_['IP9']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F100'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.855499557))
        ig = ig_['RP9']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F101'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(9.0900827198))
        ig = ig_['IP9']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F102'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-9.090082719))
        ig = ig_['RP9']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F103'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(5.3260550395))
        ig = ig_['IP9']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F104'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(24.092506375))
        ig = ig_['RP9']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F105'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-3.902049552))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F106'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(10.365394127))
        ig = ig_['IP9']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F107'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-3.902049552))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F108'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-10.36539412))
        ig = ig_['RP9']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F109'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.424005487))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F110'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(3.0290504569))
        ig = ig_['IP9']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F111'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.424005487))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F112'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-3.029050456))
        ig = ig_['RP10']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F113'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-3.902049552))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F114'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(10.365394127))
        ig = ig_['IP10']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F115'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-3.902049552))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F116'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-10.36539412))
        ig = ig_['RP10']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F117'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(5.7829343061))
        ig = ig_['IP10']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F118'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(14.768337877))
        ig = ig_['RP10']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F119'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.880884753))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F120'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(4.4029437495))
        ig = ig_['IP10']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F121'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.880884753))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F122'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-4.402943749))
        ig = ig_['RP11']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F123'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.955028563))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F124'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(4.0940743442))
        ig = ig_['IP11']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F125'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.955028563))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F126'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-4.094074344))
        ig = ig_['RP11']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F127'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.880884753))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F128'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(4.4029437495))
        ig = ig_['IP11']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F129'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.880884753))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F130'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-4.402943749))
        ig = ig_['RP11']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F131'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(3.8359133169))
        ig = ig_['IP11']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F132'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(8.4970180937))
        ig = ig_['RP12']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F133'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.525967440))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F134'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(3.175963965))
        ig = ig_['IP12']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F135'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.525967440))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F136'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-3.175963965))
        ig = ig_['RP12']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F137'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(4.0149920273))
        ig = ig_['IP12']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F138'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(5.4279385912))
        ig = ig_['RP12']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F139'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-2.489024586))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F140'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(2.2519746262))
        ig = ig_['IP12']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F141'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-2.489024586))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F142'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-2.251974626))
        ig = ig_['RP13']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F143'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-3.098927403))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F144'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(6.1027554482))
        ig = ig_['IP13']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F145'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-3.098927403))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F146'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-6.102755448))
        ig = ig_['RP13']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F147'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-2.489024586))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F148'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(2.2519746262))
        ig = ig_['IP13']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F149'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-2.489024586))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F150'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-2.251974626))
        ig = ig_['RP13']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F151'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(6.7249461485))
        ig = ig_['IP13']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F152'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(10.669693549))
        ig = ig_['RP13']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F153'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.136994157))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F154'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(2.3149634751))
        ig = ig_['IP13']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F155'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.136994157))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F156'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-2.314963475))
        ig = ig_['RP14']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F157'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.424005487))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F158'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(3.0290504569))
        ig = ig_['IP14']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F159'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.424005487))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F160'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-3.029050456))
        ig = ig_['RP14']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F161'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.136994157))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F162'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(2.3149634751))
        ig = ig_['IP14']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F163'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.136994157))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F164'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-2.314963475))
        ig = ig_['RP14']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F165'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(2.5609996448))
        ig = ig_['IP14']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F166'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(5.344013932))
        ig = ig_['FN1']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E1'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(257.14793297))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E2'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-515.1003629))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E3'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.2639541485))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E4'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(257.95312698))
        ig = ig_['FN2']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E5'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(18.779796341))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E6'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-37.76674355))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E7'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.0504741547))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E8'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(18.987552378))
        ig = ig_['FN3']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E9'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(23.945517773))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E10'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-48.09952193))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E11'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.0497138406))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E12'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(24.154483769))
        ig = ig_['FN4']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E13'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(28.840860058))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E14'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-57.85508062))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E15'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.0573251271))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E16'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(29.014509561))
        ig = ig_['FN5']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E17'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(29.691347384))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E18'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-59.56180607))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E19'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.0588594324))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E20'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(29.870757982))
        ig = ig_['FN6']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E21'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(29.572165175))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E22'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-59.20912928))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E23'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.0254204890))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E24'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(29.637005073))
        ig = ig_['FN7']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E25'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(512.43300835))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E26'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1024.866016))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E27'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(512.43300835))
        ig = ig_['FN8']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E28'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(24.995017225))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E29'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-48.89025369))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E30'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(23.907334055))
        ig = ig_['FN9']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E31'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(3.6666896805))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E32'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-7.106044600))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E33'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(3.4428786091))
        ig = ig_['FN10']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E34'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(20.86730367))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E35'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-38.89665404))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E36'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(18.125840783))
        ig = ig_['FN11']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E37'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(20.583581419))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E38'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-41.16716283))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E39'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(20.583581419))
        ig = ig_['FN12']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E40'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(12.415323736))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E41'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-24.83064747))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E42'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(12.415323736))
        ig = ig_['FN13']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E43'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(46.846975115))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E44'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-93.69395022))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E45'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(46.846975115))
        ig = ig_['FN14']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E46'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(32.22810018))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E47'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-64.45620036))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E48'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(32.22810018))
        ig = ig_['FN15']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E49'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(82.629603852))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E50'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-165.2592077))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E51'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(82.629603852))
        ig = ig_['FN16']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E52'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(122.66738612))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E53'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-245.3347722))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E54'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(122.66738612))
        ig = ig_['FN17']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E55'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(11.202938298))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E56'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-22.40587659))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E57'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(11.202938298))
        ig = ig_['FN18']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E58'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(22.923641118))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E59'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-45.84728223))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E60'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(22.923641118))
        ig = ig_['FN19']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E61'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(11.266633111))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E62'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-22.53326622))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E63'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(11.266633111))
        ig = ig_['FN20']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E64'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(6.651811606))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E65'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-13.30362321))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E66'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(6.651811606))
        ig = ig_['TN1']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E67'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(257.14793297))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E68'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-515.1003629))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E69'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.2639541485))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E70'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(257.95312698))
        ig = ig_['TN2']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E71'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(18.779796341))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E72'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-37.76674355))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E73'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.0504741547))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E74'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(18.987552378))
        ig = ig_['TN3']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E75'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(23.945517773))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E76'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-48.09952193))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E77'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.0497138406))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E78'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(24.154483769))
        ig = ig_['TN4']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E79'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(28.840860058))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E80'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-57.85508062))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E81'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.0573251271))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E82'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(29.014509561))
        ig = ig_['TN5']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E83'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(29.691347384))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E84'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-59.56180607))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E85'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.0588594324))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E86'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(29.870757982))
        ig = ig_['TN6']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E87'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(29.572165175))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E88'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-59.20912928))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E89'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.0254204890))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E90'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(29.637005073))
        ig = ig_['TN7']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E91'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(512.43300835))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E92'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1024.866016))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E93'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(512.43300835))
        ig = ig_['TN8']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E94'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(24.995017225))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E95'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-48.89025369))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E96'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(23.907334055))
        ig = ig_['TN9']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E97'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(3.6666896805))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E98'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-7.106044600))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E99'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(3.4428786091))
        ig = ig_['TN10']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E100'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(20.86730367))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E101'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-38.89665404))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E102'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(18.125840783))
        ig = ig_['TN11']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E103'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(20.583581419))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E104'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-41.16716283))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E105'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(20.583581419))
        ig = ig_['TN12']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E106'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(12.415323736))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E107'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-24.83064747))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E108'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(12.415323736))
        ig = ig_['TN13']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E109'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(46.846975115))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E110'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-93.69395022))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E111'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(46.846975115))
        ig = ig_['TN14']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E112'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(32.22810018))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E113'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-64.45620036))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E114'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(32.22810018))
        ig = ig_['TN15']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E115'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(82.629603852))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E116'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-165.2592077))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E117'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(82.629603852))
        ig = ig_['TN16']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E118'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(122.66738612))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E119'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-245.3347722))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E120'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(122.66738612))
        ig = ig_['TN17']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E121'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(11.202938298))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E122'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-22.40587659))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E123'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(11.202938298))
        ig = ig_['TN18']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E124'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(22.923641118))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E125'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-45.84728223))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E126'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(22.923641118))
        ig = ig_['TN19']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E127'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(11.266633111))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E128'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-22.53326622))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E129'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(11.266633111))
        ig = ig_['TN20']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E130'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(6.651811606))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E131'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-13.30362321))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E132'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(6.651811606))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.cupper[np.arange(pb.nle)] = np.zeros((pb.nle,1))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "QOR2-AY-38-68"
        pbm.H = pbm.H.tocsr()
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eP2(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]**2
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0e+0*EV_[0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0e+0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eP4(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]**4
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 4.0e+0*EV_[0]**3
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 12.0e+0*EV_[0]**2
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eP11(pbm,nargout,*args):

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
                H_[0,1] = 1.0e+0
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eP31(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[1]*(EV_[0]**3)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 3.0e+0*EV_[1]*EV_[0]**2
            g_[1] = EV_[0]**3
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 6.0e+0*EV_[1]*EV_[0]
                H_[0,1] = 3.0e+0*EV_[0]**2
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eP22(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = (EV_[0]*EV_[1])**2
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0e+0*EV_[0]*EV_[1]**2
            g_[1] = 2.0e+0*EV_[1]*EV_[0]**2
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 2.0e+0*EV_[1]**2
                H_[0,1] = 4.0e+0*EV_[0]*EV_[1]
                H_[1,0] = H_[0,1]
                H_[1,1] = 2.0e+0*EV_[0]**2
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eP211(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = (EV_[2]*EV_[1])*(EV_[0]**2)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0e+0*EV_[2]*EV_[1]*EV_[0]
            g_[1] = EV_[2]*EV_[0]**2
            g_[2] = EV_[1]*EV_[0]**2
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = 2.0e+0*EV_[2]*EV_[1]
                H_[0,1] = 2.0e+0*EV_[2]*EV_[0]
                H_[1,0] = H_[0,1]
                H_[0,2] = 2.0e+0*EV_[1]*EV_[0]
                H_[2,0] = H_[0,2]
                H_[1,2] = EV_[0]**2
                H_[2,1] = H_[1,2]
                H_[2,2] = 0.0e+0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eSIN11(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((3,4))
        IV_ = np.zeros(3)
        U_[0,0] = U_[0,0]+1
        U_[1,1] = U_[1,1]+1
        U_[2,2] = U_[2,2]+1
        U_[2,3] = U_[2,3]-1
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        IV_[2] = U_[2:3,:].dot(EV_)
        SINA = np.sin(IV_[2])
        COSA = np.cos(IV_[2])
        f_   = SINA*IV_[0]*IV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = SINA*IV_[1]
            g_[1] = SINA*IV_[0]
            g_[2] = COSA*IV_[0]*IV_[1]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = SINA
                H_[1,0] = H_[0,1]
                H_[2,0] = COSA*IV_[1]
                H_[0,2] = H_[2,0]
                H_[2,1] = COSA*IV_[0]
                H_[1,2] = H_[2,1]
                H_[2,2] = -SINA*IV_[0]*IV_[1]
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eCOS11(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((3,4))
        IV_ = np.zeros(3)
        U_[0,0] = U_[0,0]+1
        U_[1,1] = U_[1,1]+1
        U_[2,2] = U_[2,2]+1
        U_[2,3] = U_[2,3]-1
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        IV_[2] = U_[2:3,:].dot(EV_)
        SINA = np.sin(IV_[2])
        COSA = np.cos(IV_[2])
        f_   = COSA*IV_[0]*IV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = COSA*IV_[1]
            g_[1] = COSA*IV_[0]
            g_[2] = -SINA*IV_[0]*IV_[1]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = COSA
                H_[1,0] = H_[0,1]
                H_[2,0] = -SINA*IV_[1]
                H_[0,2] = H_[2,0]
                H_[2,1] = -SINA*IV_[0]
                H_[1,2] = H_[2,1]
                H_[2,2] = -COSA*IV_[0]*IV_[1]
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eSIN211(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((4,5))
        IV_ = np.zeros(4)
        U_[0,0] = U_[0,0]+1
        U_[1,1] = U_[1,1]+1
        U_[2,2] = U_[2,2]+1
        U_[3,3] = U_[3,3]+1
        U_[3,4] = U_[3,4]-1
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        IV_[2] = U_[2:3,:].dot(EV_)
        IV_[3] = U_[3:4,:].dot(EV_)
        SINA = np.sin(IV_[3])
        COSA = np.cos(IV_[3])
        f_   = (SINA*IV_[2]*IV_[1])*(IV_[0]**2)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0e+0*SINA*IV_[2]*IV_[1]*IV_[0]
            g_[1] = SINA*IV_[2]*IV_[0]**2
            g_[2] = SINA*IV_[1]*IV_[0]**2
            g_[3] = (COSA*IV_[2]*IV_[1])*(IV_[0]**2)
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,0] = 2.0e+0*SINA*IV_[2]*IV_[1]
                H_[0,1] = 2.0e+0*SINA*IV_[2]*IV_[0]
                H_[1,0] = H_[0,1]
                H_[0,2] = 2.0e+0*SINA*IV_[1]*IV_[0]
                H_[2,0] = H_[0,2]
                H_[1,2] = SINA*IV_[0]**2
                H_[2,1] = H_[1,2]
                H_[2,2] = 0.0e+0
                H_[3,0] = 2.0e+0*COSA*IV_[2]*IV_[1]*IV_[0]
                H_[0,3] = H_[3,0]
                H_[3,1] = (COSA*IV_[2])*(IV_[0]**2)
                H_[1,3] = H_[3,1]
                H_[3,2] = (COSA*IV_[1])*(IV_[0]**2)
                H_[2,3] = H_[3,2]
                H_[3,3] = (-SINA*IV_[2]*IV_[1])*(IV_[0]**2)
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eCOS211(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((4,5))
        IV_ = np.zeros(4)
        U_[0,0] = U_[0,0]+1
        U_[1,1] = U_[1,1]+1
        U_[2,2] = U_[2,2]+1
        U_[3,3] = U_[3,3]+1
        U_[3,4] = U_[3,4]-1
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        IV_[2] = U_[2:3,:].dot(EV_)
        IV_[3] = U_[3:4,:].dot(EV_)
        SINA = np.sin(IV_[3])
        COSA = np.cos(IV_[3])
        f_   = (COSA*IV_[2]*IV_[1])*(IV_[0]**2)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0e+0*COSA*IV_[2]*IV_[1]*IV_[0]
            g_[1] = COSA*IV_[2]*IV_[0]**2
            g_[2] = COSA*IV_[1]*IV_[0]**2
            g_[3] = (-SINA*IV_[2]*IV_[1])*(IV_[0]**2)
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,0] = 2.0e+0*COSA*IV_[2]*IV_[1]
                H_[0,1] = 2.0e+0*COSA*IV_[2]*IV_[0]
                H_[1,0] = H_[0,1]
                H_[0,2] = 2.0e+0*COSA*IV_[1]*IV_[0]
                H_[2,0] = H_[0,2]
                H_[1,2] = COSA*IV_[0]**2
                H_[2,1] = H_[1,2]
                H_[2,2] = 0.0e+0
                H_[3,0] = -2.0e+0*SINA*IV_[2]*IV_[1]*IV_[0]
                H_[0,3] = H_[3,0]
                H_[3,1] = (-SINA*IV_[2])*(IV_[0]**2)
                H_[1,3] = H_[3,1]
                H_[3,2] = (-SINA*IV_[1])*(IV_[0]**2)
                H_[2,3] = H_[3,2]
                H_[3,3] = (-COSA*IV_[2]*IV_[1])*(IV_[0]**2)
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eSIN31(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((3,4))
        IV_ = np.zeros(3)
        U_[0,0] = U_[0,0]+1
        U_[1,1] = U_[1,1]+1
        U_[2,2] = U_[2,2]+1
        U_[2,3] = U_[2,3]-1
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        IV_[2] = U_[2:3,:].dot(EV_)
        SINA = np.sin(IV_[2])
        COSA = np.cos(IV_[2])
        f_   = SINA*IV_[1]*(IV_[0]**3)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 3.0e+0*SINA*IV_[1]*IV_[0]**2
            g_[1] = SINA*IV_[0]**3
            g_[2] = COSA*IV_[1]*(IV_[0]**3)
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = 6.0e+0*SINA*IV_[1]*IV_[0]
                H_[0,1] = 3.0e+0*SINA*IV_[0]**2
                H_[1,0] = H_[0,1]
                H_[2,0] = 3.0e+0*COSA*IV_[1]*(IV_[0]**2)
                H_[0,2] = H_[2,0]
                H_[2,1] = COSA*(IV_[0]**3)
                H_[1,2] = H_[2,1]
                H_[2,2] = -SINA*IV_[1]*(IV_[0]**3)
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eCOS31(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((3,4))
        IV_ = np.zeros(3)
        U_[0,0] = U_[0,0]+1
        U_[1,1] = U_[1,1]+1
        U_[2,2] = U_[2,2]+1
        U_[2,3] = U_[2,3]-1
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        IV_[2] = U_[2:3,:].dot(EV_)
        SINA = np.sin(IV_[2])
        COSA = np.cos(IV_[2])
        f_   = COSA*IV_[1]*(IV_[0]**3)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 3.0e+0*COSA*IV_[1]*IV_[0]**2
            g_[1] = COSA*IV_[0]**3
            g_[2] = -SINA*IV_[1]*(IV_[0]**3)
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = 6.0e+0*COSA*IV_[1]*IV_[0]
                H_[0,1] = 3.0e+0*COSA*IV_[0]**2
                H_[1,0] = H_[0,1]
                H_[2,0] = -3.0e+0*SINA*IV_[1]*(IV_[0]**2)
                H_[0,2] = H_[2,0]
                H_[2,1] = -SINA*(IV_[0]**3)
                H_[1,2] = H_[2,1]
                H_[2,2] = -COSA*IV_[1]*(IV_[0]**3)
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

