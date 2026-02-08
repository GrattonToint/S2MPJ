from s2mpjlib import *
class  ACOPP30(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : ACOPP30
#    *********
# 
#    An AC Optimal Power Flow (OPF) problem for the IEEE 30 Bus
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
#    classification = "C-CQOR2-AY-72-142"
# 
#    number of nodes
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'ACOPP30'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['NODES'] = 30
        v_['LIMITS'] = 6
        v_['LINES'] = 41
        v_['1'] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['NODES'])+1):
            [iv,ix_,_] = s2mpj_ii('A'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'A'+str(I))
            [iv,ix_,_] = s2mpj_ii('M'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'M'+str(I))
        for I in range(int(v_['1']),int(v_['LIMITS'])+1):
            [iv,ix_,_] = s2mpj_ii('P'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'P'+str(I))
            [iv,ix_,_] = s2mpj_ii('Q'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'Q'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['P1']])
        valA = np.append(valA,float(200.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['P2']])
        valA = np.append(valA,float(175.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['P3']])
        valA = np.append(valA,float(300.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['P4']])
        valA = np.append(valA,float(100.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['P5']])
        valA = np.append(valA,float(300.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['P6']])
        valA = np.append(valA,float(325.0))
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
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['P1']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('IP1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'IP1')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Q1']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('RP2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'RP2')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['P2']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('IP2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'IP2')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Q2']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('RP13',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'RP13')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['P3']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('IP13',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'IP13')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Q3']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('RP22',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'RP22')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['P4']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('IP22',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'IP22')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Q4']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('RP23',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'RP23')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['P5']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('IP23',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'IP23')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Q5']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('RP27',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'RP27')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['P6']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('IP27',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'IP27')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Q6']])
        valA = np.append(valA,float(-1.0))
        for I in range(int(v_['1']),int(v_['LINES'])+1):
            [ig,ig_,_] = s2mpj_ii('FN'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'FN'+str(I))
        for I in range(int(v_['1']),int(v_['LINES'])+1):
            [ig,ig_,_] = s2mpj_ii('TN'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'TN'+str(I))
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
        self.gconst = arrset(self.gconst,ig_['RP2'],float(-0.217))
        self.gconst = arrset(self.gconst,ig_['RP3'],float(-0.024))
        self.gconst = arrset(self.gconst,ig_['RP4'],float(-0.076))
        self.gconst = arrset(self.gconst,ig_['RP7'],float(-0.228))
        self.gconst = arrset(self.gconst,ig_['RP8'],float(-0.3))
        self.gconst = arrset(self.gconst,ig_['RP10'],float(-0.058))
        self.gconst = arrset(self.gconst,ig_['RP12'],float(-0.112))
        self.gconst = arrset(self.gconst,ig_['RP14'],float(-0.062))
        self.gconst = arrset(self.gconst,ig_['RP15'],float(-0.082))
        self.gconst = arrset(self.gconst,ig_['RP16'],float(-0.035))
        self.gconst = arrset(self.gconst,ig_['RP17'],float(-0.09))
        self.gconst = arrset(self.gconst,ig_['RP18'],float(-0.032))
        self.gconst = arrset(self.gconst,ig_['RP19'],float(-0.095))
        self.gconst = arrset(self.gconst,ig_['RP20'],float(-0.022))
        self.gconst = arrset(self.gconst,ig_['RP21'],float(-0.175))
        self.gconst = arrset(self.gconst,ig_['RP23'],float(-0.032))
        self.gconst = arrset(self.gconst,ig_['RP24'],float(-0.087))
        self.gconst = arrset(self.gconst,ig_['RP26'],float(-0.035))
        self.gconst = arrset(self.gconst,ig_['RP29'],float(-0.024))
        self.gconst = arrset(self.gconst,ig_['RP30'],float(-0.106))
        self.gconst = arrset(self.gconst,ig_['IP2'],float(-0.127))
        self.gconst = arrset(self.gconst,ig_['IP3'],float(-0.012))
        self.gconst = arrset(self.gconst,ig_['IP4'],float(-0.016))
        self.gconst = arrset(self.gconst,ig_['IP7'],float(-0.109))
        self.gconst = arrset(self.gconst,ig_['IP8'],float(-0.3))
        self.gconst = arrset(self.gconst,ig_['IP10'],float(-0.02))
        self.gconst = arrset(self.gconst,ig_['IP12'],float(-0.075))
        self.gconst = arrset(self.gconst,ig_['IP14'],float(-0.016))
        self.gconst = arrset(self.gconst,ig_['IP15'],float(-0.025))
        self.gconst = arrset(self.gconst,ig_['IP16'],float(-0.018))
        self.gconst = arrset(self.gconst,ig_['IP17'],float(-0.058))
        self.gconst = arrset(self.gconst,ig_['IP18'],float(-0.009))
        self.gconst = arrset(self.gconst,ig_['IP19'],float(-0.034))
        self.gconst = arrset(self.gconst,ig_['IP20'],float(-0.007))
        self.gconst = arrset(self.gconst,ig_['IP21'],float(-0.112))
        self.gconst = arrset(self.gconst,ig_['IP23'],float(-0.016))
        self.gconst = arrset(self.gconst,ig_['IP24'],float(-0.067))
        self.gconst = arrset(self.gconst,ig_['IP26'],float(-0.023))
        self.gconst = arrset(self.gconst,ig_['IP29'],float(-0.009))
        self.gconst = arrset(self.gconst,ig_['IP30'],float(-0.019))
        self.gconst = arrset(self.gconst,ig_['FN1'],float(1.69))
        self.gconst = arrset(self.gconst,ig_['FN2'],float(1.69))
        self.gconst = arrset(self.gconst,ig_['FN3'],float(0.4225))
        self.gconst = arrset(self.gconst,ig_['FN4'],float(1.69))
        self.gconst = arrset(self.gconst,ig_['FN5'],float(1.69))
        self.gconst = arrset(self.gconst,ig_['FN6'],float(0.4225))
        self.gconst = arrset(self.gconst,ig_['FN7'],float(0.81))
        self.gconst = arrset(self.gconst,ig_['FN8'],float(0.49))
        self.gconst = arrset(self.gconst,ig_['FN9'],float(1.69))
        self.gconst = arrset(self.gconst,ig_['FN10'],float(0.1024))
        self.gconst = arrset(self.gconst,ig_['FN11'],float(0.4225))
        self.gconst = arrset(self.gconst,ig_['FN12'],float(0.1024))
        self.gconst = arrset(self.gconst,ig_['FN13'],float(0.4225))
        self.gconst = arrset(self.gconst,ig_['FN14'],float(0.4225))
        self.gconst = arrset(self.gconst,ig_['FN15'],float(0.4225))
        self.gconst = arrset(self.gconst,ig_['FN16'],float(0.4225))
        self.gconst = arrset(self.gconst,ig_['FN17'],float(0.1024))
        self.gconst = arrset(self.gconst,ig_['FN18'],float(0.1024))
        self.gconst = arrset(self.gconst,ig_['FN19'],float(0.1024))
        self.gconst = arrset(self.gconst,ig_['FN20'],float(0.0256))
        self.gconst = arrset(self.gconst,ig_['FN21'],float(0.0256))
        self.gconst = arrset(self.gconst,ig_['FN22'],float(0.0256))
        self.gconst = arrset(self.gconst,ig_['FN23'],float(0.0256))
        self.gconst = arrset(self.gconst,ig_['FN24'],float(0.1024))
        self.gconst = arrset(self.gconst,ig_['FN25'],float(0.1024))
        self.gconst = arrset(self.gconst,ig_['FN26'],float(0.1024))
        self.gconst = arrset(self.gconst,ig_['FN27'],float(0.1024))
        self.gconst = arrset(self.gconst,ig_['FN28'],float(0.1024))
        self.gconst = arrset(self.gconst,ig_['FN29'],float(0.1024))
        self.gconst = arrset(self.gconst,ig_['FN30'],float(0.0256))
        self.gconst = arrset(self.gconst,ig_['FN31'],float(0.0256))
        self.gconst = arrset(self.gconst,ig_['FN32'],float(0.0256))
        self.gconst = arrset(self.gconst,ig_['FN33'],float(0.0256))
        self.gconst = arrset(self.gconst,ig_['FN34'],float(0.0256))
        self.gconst = arrset(self.gconst,ig_['FN35'],float(0.0256))
        self.gconst = arrset(self.gconst,ig_['FN36'],float(0.4225))
        self.gconst = arrset(self.gconst,ig_['FN37'],float(0.0256))
        self.gconst = arrset(self.gconst,ig_['FN38'],float(0.0256))
        self.gconst = arrset(self.gconst,ig_['FN39'],float(0.0256))
        self.gconst = arrset(self.gconst,ig_['FN40'],float(0.1024))
        self.gconst = arrset(self.gconst,ig_['FN41'],float(0.1024))
        self.gconst = arrset(self.gconst,ig_['TN1'],float(1.69))
        self.gconst = arrset(self.gconst,ig_['TN2'],float(1.69))
        self.gconst = arrset(self.gconst,ig_['TN3'],float(0.4225))
        self.gconst = arrset(self.gconst,ig_['TN4'],float(1.69))
        self.gconst = arrset(self.gconst,ig_['TN5'],float(1.69))
        self.gconst = arrset(self.gconst,ig_['TN6'],float(0.4225))
        self.gconst = arrset(self.gconst,ig_['TN7'],float(0.81))
        self.gconst = arrset(self.gconst,ig_['TN8'],float(0.49))
        self.gconst = arrset(self.gconst,ig_['TN9'],float(1.69))
        self.gconst = arrset(self.gconst,ig_['TN10'],float(0.1024))
        self.gconst = arrset(self.gconst,ig_['TN11'],float(0.4225))
        self.gconst = arrset(self.gconst,ig_['TN12'],float(0.1024))
        self.gconst = arrset(self.gconst,ig_['TN13'],float(0.4225))
        self.gconst = arrset(self.gconst,ig_['TN14'],float(0.4225))
        self.gconst = arrset(self.gconst,ig_['TN15'],float(0.4225))
        self.gconst = arrset(self.gconst,ig_['TN16'],float(0.4225))
        self.gconst = arrset(self.gconst,ig_['TN17'],float(0.1024))
        self.gconst = arrset(self.gconst,ig_['TN18'],float(0.1024))
        self.gconst = arrset(self.gconst,ig_['TN19'],float(0.1024))
        self.gconst = arrset(self.gconst,ig_['TN20'],float(0.0256))
        self.gconst = arrset(self.gconst,ig_['TN21'],float(0.0256))
        self.gconst = arrset(self.gconst,ig_['TN22'],float(0.0256))
        self.gconst = arrset(self.gconst,ig_['TN23'],float(0.0256))
        self.gconst = arrset(self.gconst,ig_['TN24'],float(0.1024))
        self.gconst = arrset(self.gconst,ig_['TN25'],float(0.1024))
        self.gconst = arrset(self.gconst,ig_['TN26'],float(0.1024))
        self.gconst = arrset(self.gconst,ig_['TN27'],float(0.1024))
        self.gconst = arrset(self.gconst,ig_['TN28'],float(0.1024))
        self.gconst = arrset(self.gconst,ig_['TN29'],float(0.1024))
        self.gconst = arrset(self.gconst,ig_['TN30'],float(0.0256))
        self.gconst = arrset(self.gconst,ig_['TN31'],float(0.0256))
        self.gconst = arrset(self.gconst,ig_['TN32'],float(0.0256))
        self.gconst = arrset(self.gconst,ig_['TN33'],float(0.0256))
        self.gconst = arrset(self.gconst,ig_['TN34'],float(0.0256))
        self.gconst = arrset(self.gconst,ig_['TN35'],float(0.0256))
        self.gconst = arrset(self.gconst,ig_['TN36'],float(0.4225))
        self.gconst = arrset(self.gconst,ig_['TN37'],float(0.0256))
        self.gconst = arrset(self.gconst,ig_['TN38'],float(0.0256))
        self.gconst = arrset(self.gconst,ig_['TN39'],float(0.0256))
        self.gconst = arrset(self.gconst,ig_['TN40'],float(0.1024))
        self.gconst = arrset(self.gconst,ig_['TN41'],float(0.1024))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-1.0e+30)
        self.xupper = np.full((self.n,1),1.0e+30)
        self.xlower[ix_['M1']] = 0.95
        self.xupper[ix_['M1']] = 1.05
        self.xlower[ix_['M2']] = 0.95
        self.xupper[ix_['M2']] = 1.1
        self.xlower[ix_['M3']] = 0.95
        self.xupper[ix_['M3']] = 1.05
        self.xlower[ix_['M4']] = 0.95
        self.xupper[ix_['M4']] = 1.05
        self.xlower[ix_['M5']] = 0.95
        self.xupper[ix_['M5']] = 1.05
        self.xlower[ix_['M6']] = 0.95
        self.xupper[ix_['M6']] = 1.05
        self.xlower[ix_['M7']] = 0.95
        self.xupper[ix_['M7']] = 1.05
        self.xlower[ix_['M8']] = 0.95
        self.xupper[ix_['M8']] = 1.05
        self.xlower[ix_['M9']] = 0.95
        self.xupper[ix_['M9']] = 1.05
        self.xlower[ix_['M10']] = 0.95
        self.xupper[ix_['M10']] = 1.05
        self.xlower[ix_['M11']] = 0.95
        self.xupper[ix_['M11']] = 1.05
        self.xlower[ix_['M12']] = 0.95
        self.xupper[ix_['M12']] = 1.05
        self.xlower[ix_['M13']] = 0.95
        self.xupper[ix_['M13']] = 1.1
        self.xlower[ix_['M14']] = 0.95
        self.xupper[ix_['M14']] = 1.05
        self.xlower[ix_['M15']] = 0.95
        self.xupper[ix_['M15']] = 1.05
        self.xlower[ix_['M16']] = 0.95
        self.xupper[ix_['M16']] = 1.05
        self.xlower[ix_['M17']] = 0.95
        self.xupper[ix_['M17']] = 1.05
        self.xlower[ix_['M18']] = 0.95
        self.xupper[ix_['M18']] = 1.05
        self.xlower[ix_['M19']] = 0.95
        self.xupper[ix_['M19']] = 1.05
        self.xlower[ix_['M20']] = 0.95
        self.xupper[ix_['M20']] = 1.05
        self.xlower[ix_['M21']] = 0.95
        self.xupper[ix_['M21']] = 1.05
        self.xlower[ix_['M22']] = 0.95
        self.xupper[ix_['M22']] = 1.1
        self.xlower[ix_['M23']] = 0.95
        self.xupper[ix_['M23']] = 1.1
        self.xlower[ix_['M24']] = 0.95
        self.xupper[ix_['M24']] = 1.05
        self.xlower[ix_['M25']] = 0.95
        self.xupper[ix_['M25']] = 1.05
        self.xlower[ix_['M26']] = 0.95
        self.xupper[ix_['M26']] = 1.05
        self.xlower[ix_['M27']] = 0.95
        self.xupper[ix_['M27']] = 1.1
        self.xlower[ix_['M28']] = 0.95
        self.xupper[ix_['M28']] = 1.05
        self.xlower[ix_['M29']] = 0.95
        self.xupper[ix_['M29']] = 1.05
        self.xlower[ix_['M30']] = 0.95
        self.xupper[ix_['M30']] = 1.05
        self.xlower[ix_['P1']] = 0.0
        self.xupper[ix_['P1']] = 0.8
        self.xlower[ix_['P2']] = 0.0
        self.xupper[ix_['P2']] = 0.8
        self.xlower[ix_['P3']] = 0.0
        self.xupper[ix_['P3']] = 0.4
        self.xlower[ix_['P4']] = 0.0
        self.xupper[ix_['P4']] = 0.5
        self.xlower[ix_['P5']] = 0.0
        self.xupper[ix_['P5']] = 0.3
        self.xlower[ix_['P6']] = 0.0
        self.xupper[ix_['P6']] = 0.55
        self.xlower[ix_['Q1']] = -0.2
        self.xupper[ix_['Q1']] = 1.5
        self.xlower[ix_['Q2']] = -0.2
        self.xupper[ix_['Q2']] = 0.6
        self.xlower[ix_['Q3']] = -0.15
        self.xupper[ix_['Q3']] = 0.447
        self.xlower[ix_['Q4']] = -0.15
        self.xupper[ix_['Q4']] = 0.625
        self.xlower[ix_['Q5']] = -0.1
        self.xupper[ix_['Q5']] = 0.4
        self.xlower[ix_['Q6']] = -0.15
        self.xupper[ix_['Q6']] = 0.487
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(0.0))
        self.x0[ix_['M1']] = float(1.0)
        self.x0[ix_['M2']] = float(1.0)
        self.x0[ix_['M3']] = float(1.0)
        self.x0[ix_['M4']] = float(1.0)
        self.x0[ix_['M5']] = float(1.0)
        self.x0[ix_['M6']] = float(1.0)
        self.x0[ix_['M7']] = float(1.0)
        self.x0[ix_['M8']] = float(1.0)
        self.x0[ix_['M9']] = float(1.0)
        self.x0[ix_['M10']] = float(1.0)
        self.x0[ix_['M11']] = float(1.0)
        self.x0[ix_['M12']] = float(1.0)
        self.x0[ix_['M13']] = float(1.0)
        self.x0[ix_['M14']] = float(1.0)
        self.x0[ix_['M15']] = float(1.0)
        self.x0[ix_['M16']] = float(1.0)
        self.x0[ix_['M17']] = float(1.0)
        self.x0[ix_['M18']] = float(1.0)
        self.x0[ix_['M19']] = float(1.0)
        self.x0[ix_['M20']] = float(1.0)
        self.x0[ix_['M21']] = float(1.0)
        self.x0[ix_['M22']] = float(1.0)
        self.x0[ix_['M23']] = float(1.0)
        self.x0[ix_['M24']] = float(1.0)
        self.x0[ix_['M25']] = float(1.0)
        self.x0[ix_['M26']] = float(1.0)
        self.x0[ix_['M27']] = float(1.0)
        self.x0[ix_['M28']] = float(1.0)
        self.x0[ix_['M29']] = float(1.0)
        self.x0[ix_['M30']] = float(1.0)
        self.x0[ix_['P1']] = float(0.2354)
        self.x0[ix_['P2']] = float(0.6097)
        self.x0[ix_['P3']] = float(0.37)
        self.x0[ix_['P4']] = float(0.2159)
        self.x0[ix_['P5']] = float(0.192)
        self.x0[ix_['P6']] = float(0.2691)
        #%%%%%%%%%%%%%%%%%%%% QUADRATIC %%%%%%%%%%%%%%%%%%%
        irH  = np.array([],dtype=int)
        icH  = np.array([],dtype=int)
        valH = np.array([],dtype=float)
        irH  = np.append(irH,[ix_['P1']])
        icH  = np.append(icH,[ix_['P1']])
        valH = np.append(valH,float(400.0))
        irH  = np.append(irH,[ix_['P2']])
        icH  = np.append(icH,[ix_['P2']])
        valH = np.append(valH,float(350))
        irH  = np.append(irH,[ix_['P3']])
        icH  = np.append(icH,[ix_['P3']])
        valH = np.append(valH,float(500.0))
        irH  = np.append(irH,[ix_['P4']])
        icH  = np.append(icH,[ix_['P4']])
        valH = np.append(valH,float(1250.0))
        irH  = np.append(irH,[ix_['P5']])
        icH  = np.append(icH,[ix_['P5']])
        valH = np.append(valH,float(500.0))
        irH  = np.append(irH,[ix_['P6']])
        icH  = np.append(icH,[ix_['P6']])
        valH = np.append(valH,float(166.8))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eP2', iet_)
        elftv = loaset(elftv,it,0,'V1')
        [it,iet_,_] = s2mpj_ii( 'eP4', iet_)
        elftv = loaset(elftv,it,0,'V1')
        [it,iet_,_] = s2mpj_ii( 'eP22', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
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
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        ename = 'F1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F5'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F6'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F7'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F8'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F9'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F10'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F11'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F12'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F13'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F14'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F15'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F16'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F17'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F18'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F19'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F20'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F21'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F22'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F23'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F24'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F25'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F26'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F27'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F28'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F29'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F30'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F31'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F32'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F33'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F34'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F35'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F36'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F37'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F38'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F39'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F40'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F41'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F42'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F43'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F44'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F45'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F46'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F47'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F48'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F49'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F50'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F51'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F52'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F53'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F54'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F55'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F56'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F57'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F58'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F59'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F60'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F61'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F62'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F63'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F64'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F65'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F66'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F67'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F68'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F69'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F70'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F71'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F72'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F73'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F74'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F75'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F76'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F77'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F78'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F79'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F80'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F81'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F82'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F83'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F84'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F85'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F86'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F87'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F88'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F89'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F90'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F91'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F92'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F93'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F94'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F95'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F96'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F97'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F98'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F99'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F100'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F101'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F102'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F103'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F104'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F105'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F106'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F107'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F108'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F109'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F110'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F111'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F112'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F113'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F114'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F115'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F116'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M11'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A11'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F117'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M11'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A11'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F118'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F119'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F120'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F121'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F122'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F123'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F124'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F125'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F126'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F127'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F128'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F129'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F130'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F131'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F132'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F133'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F134'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F135'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F136'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F137'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F138'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F139'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F140'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M11'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A11'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F141'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M11'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A11'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F142'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M11'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F143'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F144'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F145'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F146'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F147'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M13'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A13'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F148'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M13'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A13'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F149'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F150'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F151'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F152'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F153'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F154'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F155'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F156'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F157'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F158'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F159'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F160'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F161'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M13'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A13'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F162'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M13'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A13'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F163'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M13'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F164'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F165'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F166'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F167'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F168'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F169'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F170'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F171'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F172'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F173'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F174'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F175'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F176'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F177'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F178'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F179'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F180'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F181'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F182'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F183'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F184'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F185'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F186'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F187'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F188'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F189'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F190'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F191'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F192'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F193'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F194'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F195'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F196'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F197'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F198'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F199'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F200'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F201'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F202'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F203'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F204'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F205'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F206'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F207'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F208'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F209'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F210'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F211'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F212'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F213'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F214'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F215'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F216'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F217'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F218'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F219'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F220'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F221'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F222'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F223'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F224'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F225'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F226'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F227'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F228'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F229'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F230'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F231'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F232'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F233'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F234'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F235'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F236'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F237'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F238'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F239'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F240'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F241'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F242'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F243'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F244'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F245'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F246'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F247'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F248'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F249'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F250'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F251'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F252'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F253'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F254'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F255'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F256'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F257'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F258'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F259'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F260'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F261'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F262'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F263'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F264'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F265'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F266'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F267'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F268'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F269'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F270'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F271'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F272'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F273'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F274'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F275'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F276'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F277'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F278'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F279'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F280'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F281'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F282'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F283'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F284'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F285'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F286'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F287'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F288'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F289'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F290'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F291'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F292'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F293'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F294'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F295'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F296'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M26'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A26'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F297'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M26'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A26'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F298'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M26'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A26'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F299'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M26'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A26'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F300'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F301'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F302'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F303'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F304'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M26'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A26'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F305'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M26'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A26'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F306'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M26'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A26'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F307'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M26'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A26'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F308'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M26'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F309'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M26'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F310'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F311'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F312'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F313'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F314'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F315'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F316'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F317'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F318'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F319'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F320'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F321'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F322'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F323'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F324'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F325'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F326'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F327'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F328'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F329'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F330'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F331'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F332'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F333'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F334'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F335'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F336'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F337'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F338'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F339'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F340'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F341'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F342'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F343'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F344'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F345'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F346'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F347'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F348'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F349'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F350'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F351'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F352'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F353'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F354'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN11')
        ielftype = arrset(ielftype,ie,iet_["eSIN11"])
        vname = 'M30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F355'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS11')
        ielftype = arrset(ielftype,ie,iet_["eCOS11"])
        vname = 'M30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F356'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F357'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP2')
        ielftype = arrset(ielftype,ie,iet_["eP2"])
        vname = 'M30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN31')
        ielftype = arrset(ielftype,ie,iet_["eSIN31"])
        vname = 'M1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E5'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E6'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E7'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN31')
        ielftype = arrset(ielftype,ie,iet_["eSIN31"])
        vname = 'M1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E8'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E9'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E10'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E11'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN31')
        ielftype = arrset(ielftype,ie,iet_["eSIN31"])
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E12'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E13'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E14'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E15'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E16'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E17'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E18'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN31')
        ielftype = arrset(ielftype,ie,iet_["eSIN31"])
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E19'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E20'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E21'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E22'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN31')
        ielftype = arrset(ielftype,ie,iet_["eSIN31"])
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E23'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E24'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E25'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E26'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E27'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E28'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E29'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN31')
        ielftype = arrset(ielftype,ie,iet_["eSIN31"])
        vname = 'M5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E30'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E31'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E32'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E33'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN31')
        ielftype = arrset(ielftype,ie,iet_["eSIN31"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E34'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E35'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E36'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E37'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E38'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E39'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E40'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E41'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E42'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E43'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E44'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E45'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M11'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A11'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E46'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M11'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E47'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E48'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E49'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E50'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E51'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E52'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E53'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E54'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M13'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A13'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E55'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M13'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E56'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E57'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E58'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E59'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E60'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E61'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E62'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E63'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E64'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E65'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E66'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E67'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E68'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E69'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E70'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E71'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E72'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E73'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E74'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E75'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E76'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E77'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E78'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E79'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E80'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E81'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E82'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E83'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E84'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E85'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E86'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E87'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E88'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E89'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E90'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E91'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E92'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E93'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E94'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E95'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E96'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E97'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E98'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E99'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E100'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E101'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E102'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E103'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E104'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E105'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E106'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E107'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E108'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M26'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A26'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E109'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M26'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E110'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E111'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E112'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E113'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E114'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E115'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E116'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E117'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E118'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E119'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E120'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E121'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E122'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E123'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E124'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E125'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E126'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E127'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN31')
        ielftype = arrset(ielftype,ie,iet_["eSIN31"])
        vname = 'M8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E128'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E129'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E130'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E131'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN31')
        ielftype = arrset(ielftype,ie,iet_["eSIN31"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E132'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E133'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E134'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E135'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN31')
        ielftype = arrset(ielftype,ie,iet_["eSIN31"])
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E136'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E137'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E138'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E139'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN31')
        ielftype = arrset(ielftype,ie,iet_["eSIN31"])
        vname = 'M3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A1'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E140'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E141'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E142'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E143'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN31')
        ielftype = arrset(ielftype,ie,iet_["eSIN31"])
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E144'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E145'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E146'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A3'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E147'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E148'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E149'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E150'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN31')
        ielftype = arrset(ielftype,ie,iet_["eSIN31"])
        vname = 'M5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E151'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E152'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E153'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E154'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN31')
        ielftype = arrset(ielftype,ie,iet_["eSIN31"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A2'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E155'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E156'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E157'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E158'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E159'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E160'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E161'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN31')
        ielftype = arrset(ielftype,ie,iet_["eSIN31"])
        vname = 'M7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A5'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E162'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E163'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E164'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E165'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN31')
        ielftype = arrset(ielftype,ie,iet_["eSIN31"])
        vname = 'M7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E166'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M7'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E167'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E168'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E169'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E170'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E171'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E172'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E173'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E174'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E175'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E176'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M11'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E177'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M11'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A11'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E178'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M11'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E179'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E180'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A9'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E181'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E182'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E183'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A4'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E184'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E185'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M13'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E186'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M13'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A13'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E187'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M13'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E188'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E189'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E190'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E191'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E192'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E193'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E194'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E195'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A12'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E196'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E197'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E198'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A14'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E199'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E200'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E201'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A16'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E202'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E203'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E204'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E205'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E206'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E207'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A18'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E208'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E209'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E210'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A19'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E211'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E212'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E213'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E214'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M20'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E215'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E216'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E217'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M17'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E218'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E219'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E220'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E221'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E222'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A10'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E223'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E224'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E225'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A21'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E226'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E227'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E228'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A15'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E229'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E230'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E231'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A22'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E232'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E233'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E234'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A23'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E235'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E236'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E237'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A24'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E238'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E239'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M26'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E240'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M26'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A26'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E241'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M26'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E242'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E243'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A25'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E244'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E245'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E246'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E247'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E248'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E249'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E250'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E251'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E252'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A27'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E253'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E254'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E255'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A29'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E256'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M30'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E257'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E258'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E259'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN31')
        ielftype = arrset(ielftype,ie,iet_["eSIN31"])
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A8'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E260'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E261'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP22')
        ielftype = arrset(ielftype,ie,iet_["eP22"])
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E262'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOS31')
        ielftype = arrset(ielftype,ie,iet_["eCOS31"])
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E263'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSIN31')
        ielftype = arrset(ielftype,ie,iet_["eSIN31"])
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'M6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A6'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='A2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E264'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eP4')
        ielftype = arrset(ielftype,ie,iet_["eP4"])
        vname = 'M28'
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float(-1.0e+30),float(1.0e+30),float(0.0)))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        ig = ig_['RP1']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F1'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(6.2953367876))
        ig = ig_['IP1']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F2'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(19.897279793))
        ig = ig_['RP1']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F3'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F4'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(15.0))
        ig = ig_['IP1']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F5'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F6'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-15.0))
        ig = ig_['RP1']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F7'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.295336787))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F8'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(4.9222797927))
        ig = ig_['IP1']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F9'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.295336787))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F10'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-4.922279792))
        ig = ig_['RP2']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F11'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F12'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(15.0))
        ig = ig_['IP2']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F13'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F14'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-15.0))
        ig = ig_['RP2']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F15'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(9.6892911011))
        ig = ig_['IP2']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F16'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(29.891651584))
        ig = ig_['RP2']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F17'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.846153846))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F18'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(5.2307692308))
        ig = ig_['IP2']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F19'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.846153846))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F20'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.230769230))
        ig = ig_['RP2']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F21'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.176470588))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F22'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(4.7058823529))
        ig = ig_['IP2']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F23'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.176470588))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F24'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-4.705882352))
        ig = ig_['RP2']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F25'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.666666666))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F26'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(5.0))
        ig = ig_['IP2']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F27'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.666666666))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F28'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.0))
        ig = ig_['RP3']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F29'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.295336787))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F30'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(4.9222797927))
        ig = ig_['IP3']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F31'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.295336787))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F32'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-4.922279792))
        ig = ig_['RP3']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F33'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(7.1776897287))
        ig = ig_['IP3']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F34'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(28.441691557))
        ig = ig_['RP3']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F35'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.882352941))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F36'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(23.529411765))
        ig = ig_['IP3']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F37'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.882352941))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F38'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-23.52941176))
        ig = ig_['RP4']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F39'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.846153846))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F40'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(5.2307692308))
        ig = ig_['IP4']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F41'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.846153846))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F42'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.230769230))
        ig = ig_['RP4']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F43'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.882352941))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F44'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(23.529411765))
        ig = ig_['IP4']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F45'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.882352941))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F46'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-23.52941176))
        ig = ig_['RP4']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F47'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(13.610859729))
        ig = ig_['IP4']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F48'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(56.125746606))
        ig = ig_['RP4']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F49'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.882352941))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F50'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(23.529411765))
        ig = ig_['IP4']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F51'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.882352941))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F52'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-23.52941176))
        ig = ig_['RP4']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F53'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(3.8461538462))
        ig = ig_['IP4']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F54'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-3.846153846))
        ig = ig_['RP5']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F55'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.176470588))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F56'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(4.7058823529))
        ig = ig_['IP5']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F57'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.176470588))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F58'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-4.705882352))
        ig = ig_['RP5']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F59'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(4.1350504699))
        ig = ig_['IP5']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F60'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(11.789574069))
        ig = ig_['RP5']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F61'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-2.958579881))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F62'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(7.100591716))
        ig = ig_['IP5']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F63'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-2.958579881))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F64'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-7.100591716))
        ig = ig_['RP6']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F65'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.666666666))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F66'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(5.0))
        ig = ig_['IP6']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F67'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.666666666))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F68'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.0))
        ig = ig_['RP6']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F69'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.882352941))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F70'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(23.529411765))
        ig = ig_['IP6']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F71'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.882352941))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F72'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-23.52941176))
        ig = ig_['RP6']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F73'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(22.54096159))
        ig = ig_['IP6']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F74'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(84.545346687))
        ig = ig_['RP6']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F75'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-4.109589041))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F76'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(10.95890411))
        ig = ig_['IP6']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F77'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-4.109589041))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F78'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-10.95890411))
        ig = ig_['RP6']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F79'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.882352941))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F80'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(23.529411765))
        ig = ig_['IP6']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F81'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.882352941))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F82'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-23.52941176))
        ig = ig_['RP6']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F83'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(4.7619047619))
        ig = ig_['IP6']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F84'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-4.761904761))
        ig = ig_['RP6']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F85'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.7857142857))
        ig = ig_['IP6']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F86'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.785714285))
        ig = ig_['RP6']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F87'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F88'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(15.0))
        ig = ig_['IP6']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F89'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F90'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-15.0))
        ig = ig_['RP7']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F91'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-2.958579881))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F92'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(7.100591716))
        ig = ig_['IP7']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F93'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-2.958579881))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F94'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-7.100591716))
        ig = ig_['RP7']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F95'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-4.109589041))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F96'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(10.95890411))
        ig = ig_['IP7']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F97'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-4.109589041))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F98'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-10.95890411))
        ig = ig_['RP7']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F99'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(7.0681689228))
        ig = ig_['IP7']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F100'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(18.049495826))
        ig = ig_['RP8']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F101'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.882352941))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F102'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(23.529411765))
        ig = ig_['IP8']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F103'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.882352941))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F104'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-23.52941176))
        ig = ig_['RP8']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F105'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(7.2584997302))
        ig = ig_['IP8']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F106'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(28.106567728))
        ig = ig_['RP8']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F107'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.376146789))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F108'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(4.5871559633))
        ig = ig_['IP8']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F109'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.376146789))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F110'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-4.587155963))
        ig = ig_['RP9']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F111'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(4.7619047619))
        ig = ig_['IP9']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F112'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-4.761904761))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F113'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(18.614718615))
        ig = ig_['RP9']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F114'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(9.0909090909))
        ig = ig_['IP9']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F115'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-9.090909090))
        ig = ig_['RP9']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F116'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(4.7619047619))
        ig = ig_['IP9']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F117'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-4.761904761))
        ig = ig_['RP10']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F118'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.7857142857))
        ig = ig_['IP10']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F119'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.785714285))
        ig = ig_['RP10']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F120'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(9.0909090909))
        ig = ig_['IP10']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F121'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-9.090909090))
        ig = ig_['RP10']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F122'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(13.560885291))
        ig = ig_['IP10']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F123'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(43.401934064))
        ig = ig_['RP10']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F124'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-4.109589041))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F125'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(10.95890411))
        ig = ig_['IP10']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F126'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-4.109589041))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F127'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-10.95890411))
        ig = ig_['RP10']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F128'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.724137931))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F129'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(4.0229885057))
        ig = ig_['IP10']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F130'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.724137931))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F131'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-4.022988505))
        ig = ig_['RP10']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F132'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.172413793))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F133'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(12.068965517))
        ig = ig_['IP10']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F134'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.172413793))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F135'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-12.06896551))
        ig = ig_['RP10']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F136'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-2.554744525))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F137'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(5.4744525547))
        ig = ig_['IP10']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F138'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-2.554744525))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F139'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.474452554))
        ig = ig_['RP11']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F140'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(4.7619047619))
        ig = ig_['IP11']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F141'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-4.761904761))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F142'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(4.7619047619))
        ig = ig_['RP12']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F143'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(3.8461538462))
        ig = ig_['IP12']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F144'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-3.846153846))
        ig = ig_['RP12']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F145'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(6.5455256796))
        ig = ig_['IP12']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F146'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(24.281049607))
        ig = ig_['RP12']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F147'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(7.1428571429))
        ig = ig_['IP12']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F148'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-7.142857142))
        ig = ig_['RP12']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F149'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.463414634))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F150'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(3.1707317073))
        ig = ig_['IP12']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F151'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.463414634))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F152'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-3.170731707))
        ig = ig_['RP12']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F153'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-3.211009174))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F154'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(5.9633027523))
        ig = ig_['IP12']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F155'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-3.211009174))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F156'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.963302752))
        ig = ig_['RP12']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F157'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.871101871))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F158'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(4.158004158))
        ig = ig_['IP12']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F159'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.871101871))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F160'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-4.158004158))
        ig = ig_['RP13']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F161'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(7.1428571429))
        ig = ig_['IP13']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F162'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-7.142857142))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F163'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(7.1428571429))
        ig = ig_['RP14']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F164'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.463414634))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F165'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(3.1707317073))
        ig = ig_['IP14']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F166'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.463414634))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F167'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-3.170731707))
        ig = ig_['RP14']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F168'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(3.952102417))
        ig = ig_['IP14']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F169'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(5.4331751462))
        ig = ig_['RP14']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F170'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-2.488687782))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F171'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(2.2624434389))
        ig = ig_['IP14']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F172'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-2.488687782))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F173'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-2.262443438))
        ig = ig_['RP15']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F174'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-3.211009174))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F175'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(5.9633027523))
        ig = ig_['IP15']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F176'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-3.211009174))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F177'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.963302752))
        ig = ig_['RP15']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F178'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-2.488687782))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F179'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(2.2624434389))
        ig = ig_['IP15']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F180'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-2.488687782))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F181'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-2.262443438))
        ig = ig_['RP15']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F182'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(9.5178787753))
        ig = ig_['IP15']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F183'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(15.862109828))
        ig = ig_['RP15']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F184'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.818181818))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F185'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(3.6363636364))
        ig = ig_['IP15']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F186'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.818181818))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F187'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-3.636363636))
        ig = ig_['RP15']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F188'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-2.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F189'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(4.0))
        ig = ig_['IP15']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F190'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-2.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F191'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-4.0))
        ig = ig_['RP16']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F192'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.871101871))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F193'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(4.158004158))
        ig = ig_['IP16']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F194'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.871101871))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F195'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-4.158004158))
        ig = ig_['RP16']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F196'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(3.7534548123))
        ig = ig_['IP16']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F197'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(8.6285923933))
        ig = ig_['RP16']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F198'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.882352941))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F199'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(4.4705882353))
        ig = ig_['IP16']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F200'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.882352941))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F201'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-4.470588235))
        ig = ig_['RP17']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F202'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-4.109589041))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F203'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(10.95890411))
        ig = ig_['IP17']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F204'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-4.109589041))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F205'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-10.95890411))
        ig = ig_['RP17']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F206'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.882352941))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F207'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(4.4705882353))
        ig = ig_['IP17']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F208'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.882352941))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F209'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-4.470588235))
        ig = ig_['RP17']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F210'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(5.9919419823))
        ig = ig_['IP17']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F211'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(15.429492345))
        ig = ig_['RP18']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F212'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.818181818))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F213'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(3.6363636364))
        ig = ig_['IP18']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F214'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.818181818))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F215'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-3.636363636))
        ig = ig_['RP18']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F216'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(4.7450110865))
        ig = ig_['IP18']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F217'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(9.977827051))
        ig = ig_['RP18']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F218'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-2.926829268))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F219'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(6.3414634146))
        ig = ig_['IP18']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F220'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-2.926829268))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F221'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-6.341463414))
        ig = ig_['RP19']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F222'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-2.926829268))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F223'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(6.3414634146))
        ig = ig_['IP19']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F224'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-2.926829268))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F225'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-6.341463414))
        ig = ig_['RP19']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F226'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(8.0992430614))
        ig = ig_['IP19']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F227'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(18.410428932))
        ig = ig_['RP19']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F228'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.172413793))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F229'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(12.068965517))
        ig = ig_['IP19']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F230'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.172413793))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F231'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-12.06896551))
        ig = ig_['RP20']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F232'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.724137931))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F233'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(4.0229885057))
        ig = ig_['IP20']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F234'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.724137931))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F235'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-4.022988505))
        ig = ig_['RP20']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F236'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.172413793))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F237'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(12.068965517))
        ig = ig_['IP20']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F238'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.172413793))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F239'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-12.06896551))
        ig = ig_['RP20']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F240'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(6.8965517241))
        ig = ig_['IP20']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F241'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(16.091954023))
        ig = ig_['RP21']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F242'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.172413793))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F243'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(12.068965517))
        ig = ig_['IP21']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F244'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.172413793))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F245'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-12.06896551))
        ig = ig_['RP21']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F246'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(25.172413793))
        ig = ig_['IP21']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F247'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(52.068965517))
        ig = ig_['RP21']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F248'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-20.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F249'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(40.0))
        ig = ig_['IP21']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F250'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-20.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F251'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-40.0))
        ig = ig_['RP22']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F252'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-2.554744525))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F253'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(5.4744525547))
        ig = ig_['IP22']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F254'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-2.554744525))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F255'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.474452554))
        ig = ig_['RP22']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F256'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-20.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F257'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(40.0))
        ig = ig_['IP22']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F258'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-20.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F259'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-40.0))
        ig = ig_['RP22']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F260'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(25.11884709))
        ig = ig_['IP22']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F261'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(49.320606401))
        ig = ig_['RP22']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F262'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-2.564102564))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F263'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(3.8461538462))
        ig = ig_['IP22']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F264'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-2.564102564))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F265'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-3.846153846))
        ig = ig_['RP23']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F266'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-2.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F267'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(4.0))
        ig = ig_['IP23']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F268'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-2.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F269'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-4.0))
        ig = ig_['RP23']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F270'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(3.4476614699))
        ig = ig_['IP23']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F271'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(7.0066815145))
        ig = ig_['RP23']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F272'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.447661469))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F273'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(3.0066815145))
        ig = ig_['IP23']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F274'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.447661469))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F275'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-3.006681514))
        ig = ig_['RP24']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F276'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-2.564102564))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F277'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(3.8461538462))
        ig = ig_['IP24']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F278'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-2.564102564))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F279'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-3.846153846))
        ig = ig_['RP24']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F280'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.447661469))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F281'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(3.0066815145))
        ig = ig_['IP24']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F282'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.447661469))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F283'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-3.006681514))
        ig = ig_['RP24']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F284'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(5.3221088616))
        ig = ig_['IP24']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F285'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(9.1282974296))
        ig = ig_['RP24']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F286'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.310344827))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F287'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(2.275862069))
        ig = ig_['IP24']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F288'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.310344827))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F289'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-2.275862069))
        ig = ig_['RP25']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F290'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.310344827))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F291'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(2.275862069))
        ig = ig_['IP25']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F292'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.310344827))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F293'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-2.275862069))
        ig = ig_['RP25']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F294'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(4.475953396))
        ig = ig_['IP25']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F295'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(7.8491529293))
        ig = ig_['RP25']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F296'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.208313194))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F297'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.8366360561))
        ig = ig_['IP25']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F298'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.208313194))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F299'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.836636056))
        ig = ig_['RP25']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F300'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.957295373))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F301'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(3.7366548043))
        ig = ig_['IP25']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F302'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.957295373))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F303'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-3.736654804))
        ig = ig_['RP26']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F304'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.208313194))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F305'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.8366360561))
        ig = ig_['IP26']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F306'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.208313194))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F307'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.836636056))
        ig = ig_['RP26']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F308'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.2083131948))
        ig = ig_['IP26']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F309'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.8366360561))
        ig = ig_['RP27']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F310'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.957295373))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F311'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(3.7366548043))
        ig = ig_['IP27']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F312'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.957295373))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F313'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-3.736654804))
        ig = ig_['RP27']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F314'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(3.627984583))
        ig = ig_['IP27']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F315'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(9.4025600611))
        ig = ig_['RP27']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F316'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(2.5))
        ig = ig_['IP27']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F317'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-2.5))
        ig = ig_['RP27']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F318'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-0.978647686))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F319'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.8683274021))
        ig = ig_['IP27']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F320'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-0.978647686))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F321'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.868327402))
        ig = ig_['RP27']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F322'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-0.692041522))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F323'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.2975778547))
        ig = ig_['IP27']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F324'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-0.692041522))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F325'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.297577854))
        ig = ig_['RP28']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F326'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F327'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(15.0))
        ig = ig_['IP28']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F328'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F329'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-15.0))
        ig = ig_['RP28']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F330'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.376146789))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F331'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(4.5871559633))
        ig = ig_['IP28']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F332'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.376146789))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F333'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-4.587155963))
        ig = ig_['RP28']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F334'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(2.5))
        ig = ig_['IP28']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F335'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-2.5))
        ig = ig_['RP28']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F336'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(6.376146789))
        ig = ig_['IP28']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F337'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(22.072155963))
        ig = ig_['RP29']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F338'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-0.978647686))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F339'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.8683274021))
        ig = ig_['IP29']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F340'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-0.978647686))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F341'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.868327402))
        ig = ig_['RP29']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F342'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.9013697168))
        ig = ig_['IP29']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F343'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(3.5984312084))
        ig = ig_['RP29']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F344'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-0.922722029))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F345'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.7301038062))
        ig = ig_['IP29']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F346'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-0.922722029))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F347'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.730103806))
        ig = ig_['RP30']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F348'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-0.692041522))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F349'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.2975778547))
        ig = ig_['IP30']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F350'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-0.692041522))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F351'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.297577854))
        ig = ig_['RP30']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F352'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-0.922722029))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F353'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.7301038062))
        ig = ig_['IP30']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F354'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-0.922722029))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F355'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.730103806))
        ig = ig_['RP30']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F356'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.6147635525))
        ig = ig_['IP30']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F357'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(3.0276816609))
        ig = ig_['FN1']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E1'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(249.550225))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E2'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-499.55))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E3'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.15))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E4'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(250.0))
        ig = ig_['FN2']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E5'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(25.808390155))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E6'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-51.71502590))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E7'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.0259067357))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E8'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(25.906735751))
        ig = ig_['FN3']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E9'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(30.664715385))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E10'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-61.43384615))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E11'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.0369230769))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E12'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(30.769230769))
        ig = ig_['FN4']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E13'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(588.23529412))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E14'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1176.470588))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E15'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(588.23529412))
        ig = ig_['FN5']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E16'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(23.435394118))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E17'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-46.96470588))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E18'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.0235294117))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E19'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(23.529411765))
        ig = ig_['FN6']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E20'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(27.677877778))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E21'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-55.45555555))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E22'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.0333333333))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E23'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(27.777777778))
        ig = ig_['FN7']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E24'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(588.23529412))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E25'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1176.470588))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E26'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(588.23529412))
        ig = ig_['FN8']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E27'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(59.100616716))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E28'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-118.2721893))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E29'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.0295857988))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E30'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(59.171597633))
        ig = ig_['FN9']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E31'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(136.87673733))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E32'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-273.8630137))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E33'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.0410958904))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E34'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(136.98630137))
        ig = ig_['FN10']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E35'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(588.23529412))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E36'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1176.470588))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E37'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(588.23529412))
        ig = ig_['FN11']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E38'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(22.675736961))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E39'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-45.35147392))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E40'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(22.675736961))
        ig = ig_['FN12']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E41'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(3.1887755102))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E42'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-6.377551020))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E43'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(3.1887755102))
        ig = ig_['FN13']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E44'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(22.675736961))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E45'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-45.35147392))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E46'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(22.675736961))
        ig = ig_['FN14']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E47'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(82.644628099))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E48'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-165.2892562))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E49'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(82.644628099))
        ig = ig_['FN15']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E50'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(14.792899408))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E51'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-29.58579881))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E52'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(14.792899408))
        ig = ig_['FN16']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E53'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(51.020408163))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E54'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-102.0408163))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E55'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(51.020408163))
        ig = ig_['FN17']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E56'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(12.195121951))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E57'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-24.39024390))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E58'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(12.195121951))
        ig = ig_['FN18']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E59'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(45.871559633))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E60'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-91.74311926))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E61'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(45.871559633))
        ig = ig_['FN19']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E62'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(20.79002079))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E63'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-41.58004158))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E64'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(20.79002079))
        ig = ig_['FN20']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E65'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(11.312217195))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E66'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-22.62443438))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E67'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(11.312217195))
        ig = ig_['FN21']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E68'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(23.529411765))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E69'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-47.05882352))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E70'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(23.529411765))
        ig = ig_['FN22']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E71'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(16.52892562))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E72'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-33.05785124))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E73'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(16.52892562))
        ig = ig_['FN23']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E74'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(48.780487805))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E75'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-97.56097561))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E76'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(48.780487805))
        ig = ig_['FN24']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E77'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(172.4137931))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E78'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-344.8275862))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E79'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(172.4137931))
        ig = ig_['FN25']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E80'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(19.157088123))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E81'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-38.31417624))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E82'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(19.157088123))
        ig = ig_['FN26']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E83'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(136.98630137))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E84'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-273.9726027))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E85'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(136.98630137))
        ig = ig_['FN27']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E86'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(172.4137931))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E87'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-344.8275862))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E88'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(172.4137931))
        ig = ig_['FN28']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E89'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(36.496350365))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E90'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-72.99270073))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E91'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(36.496350365))
        ig = ig_['FN29']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E92'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(2000.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E93'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-4000.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E94'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(2000.0))
        ig = ig_['FN30']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E95'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(20.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E96'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-40.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E97'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(20.0))
        ig = ig_['FN31']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E98'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(21.367521368))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E99'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-42.73504273))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E100'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(21.367521368))
        ig = ig_['FN32']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E101'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(11.135857461))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E102'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-22.27171492))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E103'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(11.135857461))
        ig = ig_['FN33']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E104'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(6.8965517241))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E105'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-13.79310344))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E106'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(6.8965517241))
        ig = ig_['FN34']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E107'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(4.8332527791))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E108'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-9.666505558))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E109'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(4.8332527791))
        ig = ig_['FN35']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E110'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(17.793594306))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E111'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-35.58718861))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E112'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(17.793594306))
        ig = ig_['FN36']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E113'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(6.25))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E114'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-12.5))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E115'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(6.25))
        ig = ig_['FN37']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E116'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(4.4483985765))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E117'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-8.896797153))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E118'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(4.4483985765))
        ig = ig_['FN38']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E119'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(2.1626297578))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E120'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-4.325259515))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E121'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(2.1626297578))
        ig = ig_['FN39']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E122'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(3.844675125))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E123'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-7.689350249))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E124'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(3.844675125))
        ig = ig_['FN40']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E125'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(22.844136697))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E126'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-45.77981651))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E127'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.0275229357))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E128'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(22.935779817))
        ig = ig_['FN41']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E129'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(249.850025))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E130'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-499.85))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E131'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.05))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E132'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(250.0))
        ig = ig_['TN1']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E133'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(249.550225))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E134'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-499.55))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E135'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.15))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E136'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(250.0))
        ig = ig_['TN2']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E137'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(25.808390155))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E138'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-51.71502590))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E139'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.0259067357))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E140'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(25.906735751))
        ig = ig_['TN3']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E141'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(30.664715385))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E142'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-61.43384615))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E143'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.0369230769))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E144'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(30.769230769))
        ig = ig_['TN4']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E145'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(588.23529412))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E146'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1176.470588))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E147'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(588.23529412))
        ig = ig_['TN5']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E148'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(23.435394118))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E149'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-46.96470588))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E150'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.0235294117))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E151'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(23.529411765))
        ig = ig_['TN6']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E152'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(27.677877778))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E153'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-55.45555555))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E154'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.0333333333))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E155'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(27.777777778))
        ig = ig_['TN7']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E156'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(588.23529412))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E157'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1176.470588))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E158'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(588.23529412))
        ig = ig_['TN8']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E159'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(59.100616716))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E160'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-118.2721893))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E161'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.0295857988))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E162'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(59.171597633))
        ig = ig_['TN9']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E163'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(136.87673733))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E164'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-273.8630137))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E165'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.0410958904))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E166'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(136.98630137))
        ig = ig_['TN10']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E167'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(588.23529412))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E168'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1176.470588))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E169'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(588.23529412))
        ig = ig_['TN11']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E170'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(22.675736961))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E171'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-45.35147392))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E172'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(22.675736961))
        ig = ig_['TN12']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E173'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(3.1887755102))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E174'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-6.377551020))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E175'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(3.1887755102))
        ig = ig_['TN13']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E176'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(22.675736961))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E177'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-45.35147392))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E178'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(22.675736961))
        ig = ig_['TN14']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E179'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(82.644628099))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E180'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-165.2892562))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E181'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(82.644628099))
        ig = ig_['TN15']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E182'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(14.792899408))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E183'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-29.58579881))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E184'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(14.792899408))
        ig = ig_['TN16']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E185'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(51.020408163))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E186'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-102.0408163))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E187'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(51.020408163))
        ig = ig_['TN17']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E188'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(12.195121951))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E189'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-24.39024390))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E190'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(12.195121951))
        ig = ig_['TN18']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E191'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(45.871559633))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E192'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-91.74311926))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E193'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(45.871559633))
        ig = ig_['TN19']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E194'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(20.79002079))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E195'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-41.58004158))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E196'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(20.79002079))
        ig = ig_['TN20']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E197'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(11.312217195))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E198'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-22.62443438))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E199'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(11.312217195))
        ig = ig_['TN21']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E200'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(23.529411765))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E201'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-47.05882352))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E202'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(23.529411765))
        ig = ig_['TN22']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E203'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(16.52892562))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E204'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-33.05785124))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E205'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(16.52892562))
        ig = ig_['TN23']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E206'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(48.780487805))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E207'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-97.56097561))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E208'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(48.780487805))
        ig = ig_['TN24']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E209'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(172.4137931))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E210'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-344.8275862))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E211'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(172.4137931))
        ig = ig_['TN25']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E212'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(19.157088123))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E213'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-38.31417624))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E214'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(19.157088123))
        ig = ig_['TN26']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E215'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(136.98630137))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E216'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-273.9726027))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E217'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(136.98630137))
        ig = ig_['TN27']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E218'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(172.4137931))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E219'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-344.8275862))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E220'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(172.4137931))
        ig = ig_['TN28']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E221'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(36.496350365))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E222'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-72.99270073))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E223'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(36.496350365))
        ig = ig_['TN29']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E224'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(2000.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E225'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-4000.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E226'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(2000.0))
        ig = ig_['TN30']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E227'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(20.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E228'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-40.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E229'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(20.0))
        ig = ig_['TN31']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E230'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(21.367521368))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E231'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-42.73504273))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E232'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(21.367521368))
        ig = ig_['TN32']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E233'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(11.135857461))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E234'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-22.27171492))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E235'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(11.135857461))
        ig = ig_['TN33']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E236'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(6.8965517241))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E237'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-13.79310344))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E238'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(6.8965517241))
        ig = ig_['TN34']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E239'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(4.8332527791))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E240'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-9.666505558))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E241'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(4.8332527791))
        ig = ig_['TN35']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E242'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(17.793594306))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E243'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-35.58718861))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E244'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(17.793594306))
        ig = ig_['TN36']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E245'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(6.25))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E246'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-12.5))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E247'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(6.25))
        ig = ig_['TN37']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E248'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(4.4483985765))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E249'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-8.896797153))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E250'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(4.4483985765))
        ig = ig_['TN38']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E251'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(2.1626297578))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E252'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-4.325259515))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E253'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(2.1626297578))
        ig = ig_['TN39']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E254'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(3.844675125))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E255'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-7.689350249))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E256'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(3.844675125))
        ig = ig_['TN40']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E257'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(22.844136697))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E258'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-45.77981651))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E259'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.0275229357))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E260'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(22.935779817))
        ig = ig_['TN41']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E261'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(249.850025))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E262'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-499.85))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E263'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.05))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E264'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(250.0))
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        self.H = csr_matrix((valH,(irH,icH)),shape=(self.n,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.cupper[np.arange(self.nle)] = np.zeros((self.nle,1))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CQOR2-AY-72-142"
        self.objderlvl = 2
        self.conderlvl = [2]


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eP2(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0,0]**2
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0e+0*EV_[0,0]
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
    def eP4(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0,0]**4
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 4.0e+0*EV_[0,0]**3
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 12.0e+0*EV_[0,0]**2
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eP22(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = (EV_[0,0]*EV_[1,0])**2
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0e+0*EV_[0,0]*EV_[1,0]**2
            g_[1] = 2.0e+0*EV_[1,0]*EV_[0,0]**2
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 2.0e+0*EV_[1,0]**2
                H_[0,1] = 4.0e+0*EV_[0,0]*EV_[1,0]
                H_[1,0] = H_[0,1]
                H_[1,1] = 2.0e+0*EV_[0,0]**2
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eSIN11(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((3,4))
        IV_ = np.zeros(3)
        U_[0,0] = U_[0,0]+1
        U_[1,1] = U_[1,1]+1
        U_[2,2] = U_[2,2]+1
        U_[2,3] = U_[2,3]-1
        IV_[0] = to_scalar(U_[0:1,:].dot(EV_))
        IV_[1] = to_scalar(U_[1:2,:].dot(EV_))
        IV_[2] = to_scalar(U_[2:3,:].dot(EV_))
        SINA = np.sin(IV_[2])
        COSA = np.cos(IV_[2])
        f_   = SINA*IV_[0]*IV_[1]
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
    def eCOS11(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((3,4))
        IV_ = np.zeros(3)
        U_[0,0] = U_[0,0]+1
        U_[1,1] = U_[1,1]+1
        U_[2,2] = U_[2,2]+1
        U_[2,3] = U_[2,3]-1
        IV_[0] = to_scalar(U_[0:1,:].dot(EV_))
        IV_[1] = to_scalar(U_[1:2,:].dot(EV_))
        IV_[2] = to_scalar(U_[2:3,:].dot(EV_))
        SINA = np.sin(IV_[2])
        COSA = np.cos(IV_[2])
        f_   = COSA*IV_[0]*IV_[1]
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
    def eSIN211(self, nargout,*args):

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
        IV_[0] = to_scalar(U_[0:1,:].dot(EV_))
        IV_[1] = to_scalar(U_[1:2,:].dot(EV_))
        IV_[2] = to_scalar(U_[2:3,:].dot(EV_))
        IV_[3] = to_scalar(U_[3:4,:].dot(EV_))
        SINA = np.sin(IV_[3])
        COSA = np.cos(IV_[3])
        f_   = (SINA*IV_[2]*IV_[1])*(IV_[0]**2)
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
    def eCOS211(self, nargout,*args):

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
        IV_[0] = to_scalar(U_[0:1,:].dot(EV_))
        IV_[1] = to_scalar(U_[1:2,:].dot(EV_))
        IV_[2] = to_scalar(U_[2:3,:].dot(EV_))
        IV_[3] = to_scalar(U_[3:4,:].dot(EV_))
        SINA = np.sin(IV_[3])
        COSA = np.cos(IV_[3])
        f_   = (COSA*IV_[2]*IV_[1])*(IV_[0]**2)
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
    def eSIN31(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((3,4))
        IV_ = np.zeros(3)
        U_[0,0] = U_[0,0]+1
        U_[1,1] = U_[1,1]+1
        U_[2,2] = U_[2,2]+1
        U_[2,3] = U_[2,3]-1
        IV_[0] = to_scalar(U_[0:1,:].dot(EV_))
        IV_[1] = to_scalar(U_[1:2,:].dot(EV_))
        IV_[2] = to_scalar(U_[2:3,:].dot(EV_))
        SINA = np.sin(IV_[2])
        COSA = np.cos(IV_[2])
        f_   = SINA*IV_[1]*(IV_[0]**3)
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
    def eCOS31(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((3,4))
        IV_ = np.zeros(3)
        U_[0,0] = U_[0,0]+1
        U_[1,1] = U_[1,1]+1
        U_[2,2] = U_[2,2]+1
        U_[2,3] = U_[2,3]-1
        IV_[0] = to_scalar(U_[0:1,:].dot(EV_))
        IV_[1] = to_scalar(U_[1:2,:].dot(EV_))
        IV_[2] = to_scalar(U_[2:3,:].dot(EV_))
        SINA = np.sin(IV_[2])
        COSA = np.cos(IV_[2])
        f_   = COSA*IV_[1]*(IV_[0]**3)
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

