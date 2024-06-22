from s2mpjlib import *
class  ANTWERP(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : ANTWERP
#    *********
# 
#    This problem arises in the determination of a synthetic population for
#    Belgian municipalities. The question is to estimate the distribution in 
#    Antwerp that households of the following types:
#       type F (a couple + 1 to 5 children + 0 to 2 additional adults)
#       type W (a woman  + 1 to 3 children + 0 to 2 additional adults)
#       type M (a man    + 1 to 3 children + 0 to 2 additional adults).
#    The data consists in 
#       - the number of individuals in households with 3 to 8 members, 
#       - the number of F, W and N households according to their number of
#         children
#       - and the total number of adults and children.
#    If we define the variables by
#    p1F: probability for a F_household to have 1 child,
#    p2F: probability for a F_household to have 2 children,
#    p3F: probability for a F_household to have 3 children,
#    p4F: probability for a F_household to have 4 children,
#    p5F: probability for a F_household to have 5 children,
#    p1W: probability for a W_household to have 1 child,
#    p2W: probability for a W_household to have 2 children,
#    p3W: probability for a W_household to have 3 children,
#    p1M: probability for a M_household to have 1 child,
#    p2M: probability for a M_household to have 2 children,
#    p3M: probability for a M_household to have 3 children,
#    q0F: probability for a F_household to have 1 additional adult,
#    q1F: probability for a F_household to have 2 additional adults,
#    q2F: probability for a F_household to have 3 additional adults,
#    q0W: probability for a W_household to have 1 additional adult,
#    q1W: probability for a W_household to have 2 additional adults,
#    q2W: probability for a W_household to have 3 additional adults,
#    q0M: probability for a M_household to have 1 additional adult,
#    q1M: probability for a M_household to have 2 additional adults,
#    q2M: probability for a M_household to have 3 additional adults,
#    nF : number of F-households,
#    nW : number of W-households,
#    nM : number of N-households,
#    nC2: number of individuals considered as children in age class 2
#    nC3: number of individuals considered as children in age class 3,
#    nA2: number of individuals considered as adults in age class 2,
#    nA3: number of individuals considered as adults in age class 3,
#    the derived predictions for the observed quantities are then given by
#    1) prediction of the number of individuals in household of size 3:
#       p1F*q0F*nF + (p1W*q1W+p2W*q0W)*nW + (p1M*q1M+p2M*q0M)*nM = M3
#    2) prediction of the number of individuals in household of size 4:
#       (p2F*q0F+p1F*q1F)*nF + (p1W*q2W+p2W*q1W+p3W*q0W)*nW 
#                            + (p1M*q2M+p2M*q1M+p3M*q0M)*nM = M4
#    3) prediction of the number of individuals in household of size 5:
#       (p3F*q0F+p2F*q1F+p1F*q2F)*nF + (p2W*q2W+p3W*q1W)*nW 
#                                    + (p2M*q2M+p3M*q1M)*nM = M5
#    4) prediction of the number of individuals in household of size 6:
#       (p4F*q0F+p3F*q1F+p2F*q2F)*nF + p3W*q2W*nW + p3M*q2M*nM = M6
#    5) prediction of the number of individuals in household of size 7:
#       (p5F*q0F+p4F*q1F+p3F*q2F)*nF = M7
#    6) prediction of the number of individuals in household of size 8:
#       (p5F*q1F+p4F*q2F)*nF = M8
#    7) prediction of the number of F-households with 1 child
#       p1F*nF*(M1F+M2F+M3F) = M1F
#    8) prediction of the number of F-households with 2 children
#       p2F*nF*(M1F+M2F+M3F) = M2F
#    9) prediction of the number of F-households with 3 children or more
#       (p3F+p4F+p5F)*nF*(M1F+M2F+M3F) = M3F
#    10) prediction of the number of W-households with 1 child
#        p1W*nW*(M1W+M2W+M3W) = M1W
#    11) prediction of the number of W-households with 2 children
#        p2W*nW*(M1W+M2W+M3W) = M2W
#    12) prediction of the number of W-households with 3 children
#        p3W*nW*(M1W+M2W+M3W) = M3W
#    13) prediction of the number of M-households with 1 child
#        p1M*nM*(M1M+M2M+M3M) = M1M
#    14) prediction of the number of M-households with 2 children
#        p2M*nM*(M1M+M2M+M3M) = M2M
#  
#    14) prediction of the number of M-households with 3 children
#        p3M*nM*(M1M+M2M+M3M) = M3M
#    15) prediction of the number of children
#        (p1F+2*p2F+3*p3F+4*P4F+5*p5F)*NF + (p1W+2*p2W+3*p3W)*nW
#                  + (p1W+2*p2W+3*p3W)*nW - nC2 -nC3 = N0 + N1
#    16) prediction of the total number of adults
#        (2*q0F+3*q1F+4*q2F)*nF + (q0W+2*q1W+3*q2W)*nW 
#                  + (q0W+2*q1W+3*q2W)*nW - nA2 -nA3 = N4
#    17) composition of age class 2
#        nC2 + nA2 = N2
#    18) composition of age class 3
#        nC3 + nA3 = N3
#    19) prediction of the total number of individuals in F-households
#        (p1F+2*p2F+3*p3F+4*P4F+5*p5F)*NF + (2*q0F+3*q1F+4*q2F)*nF = NINF
#    20) the piF are probabilities and sum up to 1
#        p1F + p2F + p3F + p4F + p5F = 1
#    21) the qiF are probabilities and sum up to 1
#        
#        q0F + q1F + q2F = 1
#    22) the piW are probabilities and sum up to 1
#        p1W + p2W + p3W = 1
#    23) the qiW are probabilities and sum up to 1
#        
#        q0W + q1W + q2W = 1
#    24) the piM are probabilities and sum up to 1
#        p1M + p2M + p3M = 1
#    25) the qiM are probabilities and sum up to 1
#        
#        q0M + q1M + q2M = 1
#    In addition, the following inequalities are imposed
#    26) the fraction of children in age class 2 exceeds that in age class 3
#         nC2/N2 >= nC3/N3
#    27) there are more adults in age class 2 than children
#         nA2 >= nC2
#    28) there are more adults in age class 3 than children
#         nA3 >= nC3
#    and the bounds on the variables are
#        0 <= piF <= 1      ( i = 1, 2, 3, 4, 5 )
#  
#        0 <= qiF <= 1      ( i = 0, 1, 2 ) 
#        0 <= piW <= 1      ( i = 1, 2, 3 )
#  
#        0 <= qiW <= 1      ( i = 0, 1, 2 ) 
#        0 <= piM <= 1      ( i = 1, 2, 3 )
#  
#        0 <= qiM <= 1      ( i = 0, 1, 2 ) 
#        nF >= 0,  nW >= 0,  nM >= 0
#        0 <= nC2 <= N2,   0 <= nA2 <= N2
#        0 <= nC3 <= N3,   0 <= nA3 <= N3
# 
#    The problem is solved as a linearly/bound  constrained nonlinear least-squares
#    problem in 27 variables.  In the least-squares formulation,
#    each equation is scaled in a proportion inverse to its right-hand side. 
# 
#    The problem appears to be very ill-conditioned.
#    SIF input: Ph. Toint, Apr 2006.
# 
#    classification = "SLR2-RN-27-8-0-3-24-0-2-0-8-0-0-0"
# 
#    Problem initial data
# 
#    Number of households according to their sizes
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'ANTWERP'

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
        v_['M3'] = 23844.0
        v_['M4'] = 16323.0
        v_['M5'] = 6613.0
        v_['M6'] = 2535.0
        v_['M7'] = 1109.0
        v_['M8'] = 1667.0
        v_['M1F'] = 16405.0
        v_['M2F'] = 13647.0
        v_['M3F'] = 9895.0
        v_['M1M'] = 4041.0
        v_['M2M'] = 1634.0
        v_['M3M'] = 637.0
        v_['M1W'] = 10966.0
        v_['M2W'] = 4566.0
        v_['M3W'] = 1921.0
        v_['N0'] = 15866.0
        v_['N1'] = 59832.0
        v_['N2'] = 61929.0
        v_['N3'] = 32321.0
        v_['N4'] = 73650.0
        v_['NINF'] = 180055.0
        v_['NINN'] = 47677.0
        v_['TMP'] = v_['M1F']+v_['M2F']
        v_['SNF'] = v_['TMP']+v_['M3F']
        v_['TMP'] = v_['M1M']+v_['M2M']
        v_['SNM'] = v_['TMP']+v_['M3M']
        v_['TMP'] = v_['M1W']+v_['M2W']
        v_['SNW'] = v_['TMP']+v_['M3W']
        v_['N2/2'] = 0.5*v_['N2']
        v_['N3/2'] = 0.5*v_['N3']
        v_['N23/2'] = v_['N2/2']+v_['N3/2']
        v_['N01'] = v_['N0']+v_['N1']
        v_['N0123'] = v_['N01']+v_['N23/2']
        v_['N234'] = v_['N4']+v_['N23/2']
        v_['SP1F'] = v_['M1F']/v_['SNF']
        v_['SP2F'] = v_['M2F']/v_['SNF']
        v_['SP1W'] = v_['M1W']/v_['SNW']
        v_['SP2W'] = v_['M2W']/v_['SNW']
        v_['SP3W'] = v_['M3W']/v_['SNW']
        v_['SP1M'] = v_['M1M']/v_['SNM']
        v_['SP2M'] = v_['M2M']/v_['SNM']
        v_['SP3M'] = v_['M3M']/v_['SNM']
        v_['1/N3'] = 1.0/v_['N3']
        v_['-1/N2'] = -1.0/v_['N2']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        pb.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2mpj_ii('P1F',ix_)
        pb.xnames=arrset(pb.xnames,iv,'P1F')
        [iv,ix_,_] = s2mpj_ii('P2F',ix_)
        pb.xnames=arrset(pb.xnames,iv,'P2F')
        [iv,ix_,_] = s2mpj_ii('P3F',ix_)
        pb.xnames=arrset(pb.xnames,iv,'P3F')
        [iv,ix_,_] = s2mpj_ii('P4F',ix_)
        pb.xnames=arrset(pb.xnames,iv,'P4F')
        [iv,ix_,_] = s2mpj_ii('P5F',ix_)
        pb.xnames=arrset(pb.xnames,iv,'P5F')
        [iv,ix_,_] = s2mpj_ii('P1W',ix_)
        pb.xnames=arrset(pb.xnames,iv,'P1W')
        [iv,ix_,_] = s2mpj_ii('P2W',ix_)
        pb.xnames=arrset(pb.xnames,iv,'P2W')
        [iv,ix_,_] = s2mpj_ii('P3W',ix_)
        pb.xnames=arrset(pb.xnames,iv,'P3W')
        [iv,ix_,_] = s2mpj_ii('P1M',ix_)
        pb.xnames=arrset(pb.xnames,iv,'P1M')
        [iv,ix_,_] = s2mpj_ii('P2M',ix_)
        pb.xnames=arrset(pb.xnames,iv,'P2M')
        [iv,ix_,_] = s2mpj_ii('P3M',ix_)
        pb.xnames=arrset(pb.xnames,iv,'P3M')
        [iv,ix_,_] = s2mpj_ii('Q0F',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q0F')
        [iv,ix_,_] = s2mpj_ii('Q1F',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q1F')
        [iv,ix_,_] = s2mpj_ii('Q2F',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q2F')
        [iv,ix_,_] = s2mpj_ii('Q0W',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q0W')
        [iv,ix_,_] = s2mpj_ii('Q1W',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q1W')
        [iv,ix_,_] = s2mpj_ii('Q2W',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q2W')
        [iv,ix_,_] = s2mpj_ii('Q0M',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q0M')
        [iv,ix_,_] = s2mpj_ii('Q1M',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q1M')
        [iv,ix_,_] = s2mpj_ii('Q2M',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q2M')
        [iv,ix_,_] = s2mpj_ii('NF',ix_)
        pb.xnames=arrset(pb.xnames,iv,'NF')
        [iv,ix_,_] = s2mpj_ii('NW',ix_)
        pb.xnames=arrset(pb.xnames,iv,'NW')
        [iv,ix_,_] = s2mpj_ii('NM',ix_)
        pb.xnames=arrset(pb.xnames,iv,'NM')
        [iv,ix_,_] = s2mpj_ii('NC2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'NC2')
        [iv,ix_,_] = s2mpj_ii('NA2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'NA2')
        [iv,ix_,_] = s2mpj_ii('NC3',ix_)
        pb.xnames=arrset(pb.xnames,iv,'NC3')
        [iv,ix_,_] = s2mpj_ii('NA3',ix_)
        pb.xnames=arrset(pb.xnames,iv,'NA3')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('HSZ3',ig_)
        gtype = arrset(gtype,ig,'<>')
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['M3']))
        [ig,ig_,_] = s2mpj_ii('HSZ4',ig_)
        gtype = arrset(gtype,ig,'<>')
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['M4']))
        [ig,ig_,_] = s2mpj_ii('HSZ5',ig_)
        gtype = arrset(gtype,ig,'<>')
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['M5']))
        [ig,ig_,_] = s2mpj_ii('HSZ6',ig_)
        gtype = arrset(gtype,ig,'<>')
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['M6']))
        [ig,ig_,_] = s2mpj_ii('HSZ7',ig_)
        gtype = arrset(gtype,ig,'<>')
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['M7']))
        [ig,ig_,_] = s2mpj_ii('HSZ8',ig_)
        gtype = arrset(gtype,ig,'<>')
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['M8']))
        [ig,ig_,_] = s2mpj_ii('HST1F',ig_)
        gtype = arrset(gtype,ig,'<>')
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['M1F']))
        [ig,ig_,_] = s2mpj_ii('HST2F',ig_)
        gtype = arrset(gtype,ig,'<>')
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['M2F']))
        [ig,ig_,_] = s2mpj_ii('HST3F',ig_)
        gtype = arrset(gtype,ig,'<>')
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['M3F']))
        [ig,ig_,_] = s2mpj_ii('HST1W',ig_)
        gtype = arrset(gtype,ig,'<>')
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['M1W']))
        [ig,ig_,_] = s2mpj_ii('HST2W',ig_)
        gtype = arrset(gtype,ig,'<>')
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['M2W']))
        [ig,ig_,_] = s2mpj_ii('HST3W',ig_)
        gtype = arrset(gtype,ig,'<>')
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['M3W']))
        [ig,ig_,_] = s2mpj_ii('HST1M',ig_)
        gtype = arrset(gtype,ig,'<>')
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['M1M']))
        [ig,ig_,_] = s2mpj_ii('HST2M',ig_)
        gtype = arrset(gtype,ig,'<>')
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['M2M']))
        [ig,ig_,_] = s2mpj_ii('HST3M',ig_)
        gtype = arrset(gtype,ig,'<>')
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['M3M']))
        v_['WHCH'] = 100.0*v_['N0123']
        [ig,ig_,_] = s2mpj_ii('HCH',ig_)
        gtype = arrset(gtype,ig,'<>')
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['WHCH']))
        iv = ix_['NC2']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['NC3']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        v_['WHAD'] = 100.0*v_['N234']
        [ig,ig_,_] = s2mpj_ii('HAD',ig_)
        gtype = arrset(gtype,ig,'<>')
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['WHAD']))
        iv = ix_['NA2']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['NA3']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('AGE2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'AGE2')
        iv = ix_['NC2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['NA2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('AGE3',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'AGE3')
        iv = ix_['NC3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['NA3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('HINF',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('HINN',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('PSF',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PSF')
        iv = ix_['P1F']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['P2F']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['P3F']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['P4F']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['P5F']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('PSW',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PSW')
        iv = ix_['P1W']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['P2W']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['P3W']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('PSM',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PSM')
        iv = ix_['P1M']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['P2M']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['P3M']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('QSF',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'QSF')
        iv = ix_['Q0F']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['Q1F']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['Q2F']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('QSM',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'QSM')
        iv = ix_['Q0M']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['Q1M']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['Q2M']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('QSW',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'QSW')
        iv = ix_['Q0W']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['Q1W']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['Q2W']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('INEQ2',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'INEQ2')
        iv = ix_['NC2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['NA2']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('INEQ3',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'INEQ3')
        iv = ix_['NC3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['NA3']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
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
        pbm.gconst = arrset(pbm.gconst,ig_['HSZ3'],float(v_['M3']))
        pbm.gconst = arrset(pbm.gconst,ig_['HSZ4'],float(v_['M4']))
        pbm.gconst = arrset(pbm.gconst,ig_['HSZ5'],float(v_['M5']))
        pbm.gconst = arrset(pbm.gconst,ig_['HSZ6'],float(v_['M6']))
        pbm.gconst = arrset(pbm.gconst,ig_['HSZ7'],float(v_['M7']))
        pbm.gconst = arrset(pbm.gconst,ig_['HSZ8'],float(v_['M8']))
        pbm.gconst = arrset(pbm.gconst,ig_['HST1F'],float(v_['M1F']))
        pbm.gconst = arrset(pbm.gconst,ig_['HST2F'],float(v_['M2F']))
        pbm.gconst = arrset(pbm.gconst,ig_['HST3F'],float(v_['M3F']))
        pbm.gconst = arrset(pbm.gconst,ig_['HST1W'],float(v_['M1W']))
        pbm.gconst = arrset(pbm.gconst,ig_['HST2W'],float(v_['M2W']))
        pbm.gconst = arrset(pbm.gconst,ig_['HST3W'],float(v_['M3W']))
        pbm.gconst = arrset(pbm.gconst,ig_['HST1M'],float(v_['M1M']))
        pbm.gconst = arrset(pbm.gconst,ig_['HST2M'],float(v_['M2M']))
        pbm.gconst = arrset(pbm.gconst,ig_['HST3M'],float(v_['M3M']))
        pbm.gconst = arrset(pbm.gconst,ig_['HCH'],float(v_['N01']))
        pbm.gconst = arrset(pbm.gconst,ig_['HAD'],float(v_['N4']))
        pbm.gconst = arrset(pbm.gconst,ig_['HINF'],float(v_['NINF']))
        pbm.gconst = arrset(pbm.gconst,ig_['HINN'],float(v_['NINN']))
        pbm.gconst = arrset(pbm.gconst,ig_['AGE2'],float(v_['N2']))
        pbm.gconst = arrset(pbm.gconst,ig_['AGE3'],float(v_['N3']))
        pbm.gconst = arrset(pbm.gconst,ig_['PSF'],float(1.0))
        pbm.gconst = arrset(pbm.gconst,ig_['PSM'],float(1.0))
        pbm.gconst = arrset(pbm.gconst,ig_['PSW'],float(1.0))
        pbm.gconst = arrset(pbm.gconst,ig_['QSF'],float(1.0))
        pbm.gconst = arrset(pbm.gconst,ig_['QSM'],float(1.0))
        pbm.gconst = arrset(pbm.gconst,ig_['QSW'],float(1.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),0.0)
        pb.xupper = np.full((pb.n,1),1.0)
        pb.xupper[ix_['NF']] = +float('Inf')
        pb.xupper[ix_['NM']] = +float('Inf')
        pb.xupper[ix_['NW']] = +float('Inf')
        pb.xlower[ix_['NC2']] = 0.0
        pb.xupper[ix_['NC2']] = v_['N2']
        pb.xlower[ix_['NA2']] = 0.0
        pb.xupper[ix_['NA2']] = v_['N2']
        pb.xlower[ix_['NC3']] = 0.0
        pb.xupper[ix_['NC3']] = v_['N3']
        pb.xlower[ix_['NA3']] = 0.0
        pb.xupper[ix_['NA3']] = v_['N3']
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        if('P1F' in ix_):
            pb.x0[ix_['P1F']] = float(v_['SP1F'])
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['P1F']),float(v_['SP1F'])))
        if('P2F' in ix_):
            pb.x0[ix_['P2F']] = float(v_['SP2F'])
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['P2F']),float(v_['SP2F'])))
        if('P3F' in ix_):
            pb.x0[ix_['P3F']] = float(0.15)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['P3F']),float(0.15)))
        if('P4F' in ix_):
            pb.x0[ix_['P4F']] = float(0.10)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['P4F']),float(0.10)))
        if('P5F' in ix_):
            pb.x0[ix_['P5F']] = float(0.05)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['P5F']),float(0.05)))
        if('P1W' in ix_):
            pb.x0[ix_['P1W']] = float(v_['SP1W'])
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['P1W']),float(v_['SP1W'])))
        if('P2W' in ix_):
            pb.x0[ix_['P2W']] = float(v_['SP2W'])
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['P2W']),float(v_['SP2W'])))
        if('P3W' in ix_):
            pb.x0[ix_['P3W']] = float(v_['SP3W'])
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['P3W']),float(v_['SP3W'])))
        if('P1M' in ix_):
            pb.x0[ix_['P1M']] = float(v_['SP1M'])
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['P1M']),float(v_['SP1M'])))
        if('P2M' in ix_):
            pb.x0[ix_['P2M']] = float(v_['SP2M'])
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['P2M']),float(v_['SP2M'])))
        if('P3M' in ix_):
            pb.x0[ix_['P3M']] = float(v_['SP3M'])
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['P3M']),float(v_['SP3M'])))
        if('Q0F' in ix_):
            pb.x0[ix_['Q0F']] = float(0.6)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['Q0F']),float(0.6)))
        if('Q1F' in ix_):
            pb.x0[ix_['Q1F']] = float(0.3)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['Q1F']),float(0.3)))
        if('Q2F' in ix_):
            pb.x0[ix_['Q2F']] = float(0.1)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['Q2F']),float(0.1)))
        if('Q0M' in ix_):
            pb.x0[ix_['Q0M']] = float(0.6)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['Q0M']),float(0.6)))
        if('Q1M' in ix_):
            pb.x0[ix_['Q1M']] = float(0.3)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['Q1M']),float(0.3)))
        if('Q2M' in ix_):
            pb.x0[ix_['Q2M']] = float(0.1)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['Q2M']),float(0.1)))
        if('Q0W' in ix_):
            pb.x0[ix_['Q0W']] = float(0.6)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['Q0W']),float(0.6)))
        if('Q1W' in ix_):
            pb.x0[ix_['Q1W']] = float(0.3)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['Q1W']),float(0.3)))
        if('Q2W' in ix_):
            pb.x0[ix_['Q2W']] = float(0.1)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['Q2W']),float(0.1)))
        if('NF' in ix_):
            pb.x0[ix_['NF']] = float(v_['SNF'])
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['NF']),float(v_['SNF'])))
        if('NW' in ix_):
            pb.x0[ix_['NW']] = float(v_['SNW'])
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['NW']),float(v_['SNW'])))
        if('NM' in ix_):
            pb.x0[ix_['NM']] = float(v_['SNM'])
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['NM']),float(v_['SNM'])))
        if('NC2' in ix_):
            pb.x0[ix_['NC2']] = float(0.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['NC2']),float(0.0)))
        if('NC3' in ix_):
            pb.x0[ix_['NC3']] = float(0.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['NC3']),float(0.0)))
        if('NA2' in ix_):
            pb.x0[ix_['NA2']] = float(v_['N2'])
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['NA2']),float(v_['N2'])))
        if('NA3' in ix_):
            pb.x0[ix_['NA3']] = float(v_['N3'])
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['NA3']),float(v_['N3'])))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        [it,iet_,_] = s2mpj_ii( 'en3PR', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        ename = 'P1FQ0FNF'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en3PR')
        ielftype = arrset(ielftype, ie, iet_["en3PR"])
        vname = 'P1F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Q0F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NF'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P1WQ1WNW'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en3PR')
        ielftype = arrset(ielftype, ie, iet_["en3PR"])
        vname = 'P1W'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Q1W'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NW'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P2WQ0WNW'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en3PR')
        ielftype = arrset(ielftype, ie, iet_["en3PR"])
        vname = 'P2W'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Q0W'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NW'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P1MQ1MNM'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en3PR')
        ielftype = arrset(ielftype, ie, iet_["en3PR"])
        vname = 'P1M'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Q1M'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NM'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P2MQ0MNM'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en3PR')
        ielftype = arrset(ielftype, ie, iet_["en3PR"])
        vname = 'P2M'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Q0M'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NM'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P2FQ0FNF'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en3PR')
        ielftype = arrset(ielftype, ie, iet_["en3PR"])
        vname = 'P2F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Q0F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NF'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P1FQ1FNF'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en3PR')
        ielftype = arrset(ielftype, ie, iet_["en3PR"])
        vname = 'P1F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Q1F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NF'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P1WQ2WNW'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en3PR')
        ielftype = arrset(ielftype, ie, iet_["en3PR"])
        vname = 'P1W'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Q2W'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NW'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P2WQ1WNW'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en3PR')
        ielftype = arrset(ielftype, ie, iet_["en3PR"])
        vname = 'P2W'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Q1W'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NW'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P3WQ0WNW'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en3PR')
        ielftype = arrset(ielftype, ie, iet_["en3PR"])
        vname = 'P3W'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Q0W'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NW'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P1MQ2MNM'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en3PR')
        ielftype = arrset(ielftype, ie, iet_["en3PR"])
        vname = 'P1M'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Q2M'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NM'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P2MQ1MNM'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en3PR')
        ielftype = arrset(ielftype, ie, iet_["en3PR"])
        vname = 'P2M'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Q1M'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NM'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P3MQ0MNM'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en3PR')
        ielftype = arrset(ielftype, ie, iet_["en3PR"])
        vname = 'P3M'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Q0M'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NM'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P3FQ0FNF'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en3PR')
        ielftype = arrset(ielftype, ie, iet_["en3PR"])
        vname = 'P3F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Q0F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NF'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P2FQ1FNF'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en3PR')
        ielftype = arrset(ielftype, ie, iet_["en3PR"])
        vname = 'P2F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Q1F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NF'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P1FQ2FNF'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en3PR')
        ielftype = arrset(ielftype, ie, iet_["en3PR"])
        vname = 'P1F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Q2F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NF'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P2WQ2WNW'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en3PR')
        ielftype = arrset(ielftype, ie, iet_["en3PR"])
        vname = 'P2W'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Q2W'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NW'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P3WQ1WNW'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en3PR')
        ielftype = arrset(ielftype, ie, iet_["en3PR"])
        vname = 'P3W'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Q1W'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NW'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P2MQ2MNM'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en3PR')
        ielftype = arrset(ielftype, ie, iet_["en3PR"])
        vname = 'P2M'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Q2M'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NM'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P3MQ1MNM'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en3PR')
        ielftype = arrset(ielftype, ie, iet_["en3PR"])
        vname = 'P3M'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Q1M'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NM'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P4FQ0FNF'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en3PR')
        ielftype = arrset(ielftype, ie, iet_["en3PR"])
        vname = 'P4F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Q0F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NF'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P3FQ1FNF'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en3PR')
        ielftype = arrset(ielftype, ie, iet_["en3PR"])
        vname = 'P3F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Q1F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NF'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P2FQ2FNF'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en3PR')
        ielftype = arrset(ielftype, ie, iet_["en3PR"])
        vname = 'P2F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Q2F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NF'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P3WQ2WNW'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en3PR')
        ielftype = arrset(ielftype, ie, iet_["en3PR"])
        vname = 'P3W'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Q2W'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NW'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P3MQ2MNM'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en3PR')
        ielftype = arrset(ielftype, ie, iet_["en3PR"])
        vname = 'P3M'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Q2M'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NM'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P5FQ0FNF'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en3PR')
        ielftype = arrset(ielftype, ie, iet_["en3PR"])
        vname = 'P5F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Q0F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NF'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P4FQ1FNF'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en3PR')
        ielftype = arrset(ielftype, ie, iet_["en3PR"])
        vname = 'P4F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Q1F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NF'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P3FQ2FNF'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en3PR')
        ielftype = arrset(ielftype, ie, iet_["en3PR"])
        vname = 'P3F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Q2F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NF'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P5FQ1FNF'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en3PR')
        ielftype = arrset(ielftype, ie, iet_["en3PR"])
        vname = 'P5F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Q1F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NF'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P4FQ2FNF'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en3PR')
        ielftype = arrset(ielftype, ie, iet_["en3PR"])
        vname = 'P4F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Q2F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NF'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P1FNF'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
        ielftype = arrset(ielftype, ie, iet_["en2PR"])
        vname = 'P1F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NF'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P2FNF'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
        ielftype = arrset(ielftype, ie, iet_["en2PR"])
        vname = 'P2F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NF'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P3FNF'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
        ielftype = arrset(ielftype, ie, iet_["en2PR"])
        vname = 'P3F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NF'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P4FNF'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
        ielftype = arrset(ielftype, ie, iet_["en2PR"])
        vname = 'P4F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NF'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P5FNF'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
        ielftype = arrset(ielftype, ie, iet_["en2PR"])
        vname = 'P5F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NF'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P1WNW'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
        ielftype = arrset(ielftype, ie, iet_["en2PR"])
        vname = 'P1W'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NW'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P2WNW'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
        ielftype = arrset(ielftype, ie, iet_["en2PR"])
        vname = 'P2W'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NW'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P3WNW'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
        ielftype = arrset(ielftype, ie, iet_["en2PR"])
        vname = 'P3W'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NW'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P1MNM'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
        ielftype = arrset(ielftype, ie, iet_["en2PR"])
        vname = 'P1M'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NM'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P2MNM'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
        ielftype = arrset(ielftype, ie, iet_["en2PR"])
        vname = 'P2M'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NM'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'P3MNM'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
        ielftype = arrset(ielftype, ie, iet_["en2PR"])
        vname = 'P3M'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NM'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'Q0FNF'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
        ielftype = arrset(ielftype, ie, iet_["en2PR"])
        vname = 'Q0F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NF'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'Q1FNF'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
        ielftype = arrset(ielftype, ie, iet_["en2PR"])
        vname = 'Q1F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NF'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'Q2FNF'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
        ielftype = arrset(ielftype, ie, iet_["en2PR"])
        vname = 'Q2F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NF'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'Q0WNW'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
        ielftype = arrset(ielftype, ie, iet_["en2PR"])
        vname = 'Q0W'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NW'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'Q1WNW'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
        ielftype = arrset(ielftype, ie, iet_["en2PR"])
        vname = 'Q1W'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NW'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'Q2WNW'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
        ielftype = arrset(ielftype, ie, iet_["en2PR"])
        vname = 'Q2W'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NW'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'Q0MNM'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
        ielftype = arrset(ielftype, ie, iet_["en2PR"])
        vname = 'Q0M'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NM'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'Q1MNM'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
        ielftype = arrset(ielftype, ie, iet_["en2PR"])
        vname = 'Q1M'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NM'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'Q2MNM'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
        ielftype = arrset(ielftype, ie, iet_["en2PR"])
        vname = 'Q2M'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NM'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gL2',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        ig = ig_['HSZ3']
        pbm.grftype = arrset(pbm.grftype,ig,'gL2')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P1FQ0FNF'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P1WQ1WNW'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P2WQ0WNW'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P1MQ1MNM'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P2MQ0MNM'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['HSZ4']
        pbm.grftype = arrset(pbm.grftype,ig,'gL2')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P2FQ0FNF'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P1FQ1FNF'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P1WQ2WNW'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P2WQ1WNW'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P3WQ0WNW'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P1MQ2MNM'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P2MQ1MNM'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P3MQ0MNM'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_['HSZ5']
        pbm.grftype = arrset(pbm.grftype,ig,'gL2')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P3FQ0FNF'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P2FQ1FNF'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P1FQ2FNF'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P2WQ2WNW'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P3WQ1WNW'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P2MQ2MNM'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P3MQ1MNM'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['HSZ6']
        pbm.grftype = arrset(pbm.grftype,ig,'gL2')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P4FQ0FNF'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P3FQ1FNF'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P2FQ2FNF'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P3WQ2WNW'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P3MQ2MNM'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['HSZ7']
        pbm.grftype = arrset(pbm.grftype,ig,'gL2')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P5FQ0FNF'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P4FQ1FNF'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P3FQ2FNF'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['HSZ8']
        pbm.grftype = arrset(pbm.grftype,ig,'gL2')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P5FQ1FNF'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P4FQ2FNF'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_['HST1F']
        pbm.grftype = arrset(pbm.grftype,ig,'gL2')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P1FNF'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['HST2F']
        pbm.grftype = arrset(pbm.grftype,ig,'gL2')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P2FNF'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['HST3F']
        pbm.grftype = arrset(pbm.grftype,ig,'gL2')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P3FNF'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P4FNF'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P5FNF'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['HST1W']
        pbm.grftype = arrset(pbm.grftype,ig,'gL2')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P1WNW'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['HST2W']
        pbm.grftype = arrset(pbm.grftype,ig,'gL2')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P2WNW'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['HST3W']
        pbm.grftype = arrset(pbm.grftype,ig,'gL2')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P3WNW'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['HST1M']
        pbm.grftype = arrset(pbm.grftype,ig,'gL2')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P1MNM'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['HST2M']
        pbm.grftype = arrset(pbm.grftype,ig,'gL2')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P2MNM'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['HST3M']
        pbm.grftype = arrset(pbm.grftype,ig,'gL2')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P3MNM'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['HCH']
        pbm.grftype = arrset(pbm.grftype,ig,'gL2')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P1FNF'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P2FNF'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(2.0))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P3FNF'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(3.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P4FNF'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(4.0))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P5FNF'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(5.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P1MNM'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P2MNM'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(2.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P3MNM'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(3.0))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P1WNW'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P2WNW'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(2.0))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P3WNW'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(3.0))
        ig = ig_['HAD']
        pbm.grftype = arrset(pbm.grftype,ig,'gL2')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Q0FNF'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(2.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Q1FNF'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(3.0))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Q2FNF'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(4.0))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Q0MNM'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Q1MNM'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(2.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Q2MNM'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(3.0))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Q0WNW'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Q1WNW'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(2.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Q2WNW'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(3.0))
        ig = ig_['HINF']
        pbm.grftype = arrset(pbm.grftype,ig,'gL2')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P1FNF'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P2FNF'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(2.0))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P3FNF'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(3.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P4FNF'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(4.0))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P5FNF'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(5.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Q0FNF'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(2.0))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Q1FNF'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(3.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Q2FNF'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(4.0))
        ig = ig_['HINN']
        pbm.grftype = arrset(pbm.grftype,ig,'gL2')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P1WNW'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P2WNW'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(2.0))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P3WNW'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(3.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Q0WNW'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Q1WNW'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(2.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Q2WNW'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(3.0))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P1MNM'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P2MNM'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(2.0))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P3MNM'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(3.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Q0MNM'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Q1MNM'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(2.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Q2MNM'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(3.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN                0.0
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
        pb.pbclass = "SLR2-RN-27-8-0-3-24-0-2-0-8-0-0-0"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def en2PR(pbm,nargout,*args):

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

    @staticmethod
    def en3PR(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[1]*EV_[2]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]*EV_[2]
            g_[1] = EV_[0]*EV_[2]
            g_[2] = EV_[0]*EV_[1]
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = EV_[2]
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1]
                H_[2,0] = H_[0,2]
                H_[1,2] = EV_[0]
                H_[2,1] = H_[1,2]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gL2(pbm,nargout,*args):

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

