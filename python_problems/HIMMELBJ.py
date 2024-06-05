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
#    classification = "OLR2-MY-45-14"
# 
#    Number of variable sets
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HIMMELBJ'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'HIMMELBJ'
        pbm.name = self.name
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
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for K in range(int(v_['1']),int(v_['NSETS'])+1):
            v_['RNSK'] = v_['NS'+str(K)]
            v_['NSK'] = int(np.fix(v_['RNSK']))
            for J in range(int(v_['1']),int(v_['NSK'])+1):
                [iv,ix_,_] = s2mpj_ii('X'+str(J)+','+str(K),ix_)
                pb.xnames=arrset(pb.xnames,iv,'X'+str(J)+','+str(K))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X2,1']
        pbm.A[ig,iv] = float(-7.69)+pbm.A[ig,iv]
        iv = ix_['X3,1']
        pbm.A[ig,iv] = float(-11.52)+pbm.A[ig,iv]
        iv = ix_['X4,1']
        pbm.A[ig,iv] = float(-36.60)+pbm.A[ig,iv]
        iv = ix_['X1,2']
        pbm.A[ig,iv] = float(-10.94)+pbm.A[ig,iv]
        iv = ix_['X8,2']
        pbm.A[ig,iv] = float(2.5966)+pbm.A[ig,iv]
        iv = ix_['X9,2']
        pbm.A[ig,iv] = float(-39.39)+pbm.A[ig,iv]
        iv = ix_['X10,2']
        pbm.A[ig,iv] = float(-21.35)+pbm.A[ig,iv]
        iv = ix_['X11,2']
        pbm.A[ig,iv] = float(-32.84)+pbm.A[ig,iv]
        iv = ix_['X12,2']
        pbm.A[ig,iv] = float(6.26)+pbm.A[ig,iv]
        iv = ix_['X1,3']
        pbm.A[ig,iv] = float(10.45)+pbm.A[ig,iv]
        iv = ix_['X3,3']
        pbm.A[ig,iv] = float(-0.5)+pbm.A[ig,iv]
        iv = ix_['X7,3']
        pbm.A[ig,iv] = float(2.2435)+pbm.A[ig,iv]
        iv = ix_['X9,3']
        pbm.A[ig,iv] = float(-39.39)+pbm.A[ig,iv]
        iv = ix_['X10,3']
        pbm.A[ig,iv] = float(-21.49)+pbm.A[ig,iv]
        iv = ix_['X11,3']
        pbm.A[ig,iv] = float(-32.84)+pbm.A[ig,iv]
        iv = ix_['X12,3']
        pbm.A[ig,iv] = float(6.12)+pbm.A[ig,iv]
        iv = ix_['X15,3']
        pbm.A[ig,iv] = float(-1.9028)+pbm.A[ig,iv]
        iv = ix_['X16,3']
        pbm.A[ig,iv] = float(-2.8889)+pbm.A[ig,iv]
        iv = ix_['X17,3']
        pbm.A[ig,iv] = float(-3.3622)+pbm.A[ig,iv]
        iv = ix_['X18,3']
        pbm.A[ig,iv] = float(-7.4854)+pbm.A[ig,iv]
        iv = ix_['X1,4']
        pbm.A[ig,iv] = float(-15.639)+pbm.A[ig,iv]
        iv = ix_['X3,4']
        pbm.A[ig,iv] = float(21.81)+pbm.A[ig,iv]
        iv = ix_['X1,5']
        pbm.A[ig,iv] = float(-16.79)+pbm.A[ig,iv]
        iv = ix_['X3,5']
        pbm.A[ig,iv] = float(18.9779)+pbm.A[ig,iv]
        iv = ix_['X2,6']
        pbm.A[ig,iv] = float(11.959)+pbm.A[ig,iv]
        iv = ix_['X2,7']
        pbm.A[ig,iv] = float(12.899)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C1')
        iv = ix_['X1,1']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X1,2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X1,3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X15,3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X16,3']
        pbm.A[ig,iv] = float(2.0)+pbm.A[ig,iv]
        iv = ix_['X17,3']
        pbm.A[ig,iv] = float(3.0)+pbm.A[ig,iv]
        iv = ix_['X18,3']
        pbm.A[ig,iv] = float(4.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C2')
        iv = ix_['X2,1']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X2,2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X10,2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X11,2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X12,2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X2,3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X10,3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X11,3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X12,3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X2,6']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X2,7']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C3',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C3')
        iv = ix_['X3,1']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X3,2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X3,3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C4',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C4')
        iv = ix_['X4,1']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X4,2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X9,2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X11,2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X12,2']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X4,3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X9,3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X11,3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X12,3']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X1,4']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X3,4']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X1,5']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X3,5']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X2,6']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X2,7']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C5',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C5')
        iv = ix_['X4,1']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X5,2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X9,2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X10,2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X11,2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X12,2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X5,3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X9,3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X10,3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X11,3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X12,3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C6',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C6')
        iv = ix_['X6,2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X6,3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C7',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C7')
        iv = ix_['X7,2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X7,3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C8',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C8')
        iv = ix_['X8,2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X8,3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C11',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C11')
        iv = ix_['X14,3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X15,3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X16,3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X17,3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X18,3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C12',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C12')
        iv = ix_['X4,2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X5,2']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X6,2']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X7,2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X8,2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X10,2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X12,2']
        pbm.A[ig,iv] = float(-2.0)+pbm.A[ig,iv]
        iv = ix_['X13,2']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X14,3']
        pbm.A[ig,iv] = float(-4.0)+pbm.A[ig,iv]
        iv = ix_['X15,3']
        pbm.A[ig,iv] = float(-3.0)+pbm.A[ig,iv]
        iv = ix_['X16,3']
        pbm.A[ig,iv] = float(-2.0)+pbm.A[ig,iv]
        iv = ix_['X17,3']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C13',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C13')
        iv = ix_['X15,3']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X16,3']
        pbm.A[ig,iv] = float(-2.0)+pbm.A[ig,iv]
        iv = ix_['X17,3']
        pbm.A[ig,iv] = float(-3.0)+pbm.A[ig,iv]
        iv = ix_['X18,3']
        pbm.A[ig,iv] = float(-4.0)+pbm.A[ig,iv]
        iv = ix_['X1,4']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X2,4']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X3,4']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C14',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C14')
        iv = ix_['X1,5']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X2,5']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X3,5']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C15',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C15')
        iv = ix_['X3,4']
        pbm.A[ig,iv] = float(-4.0)+pbm.A[ig,iv]
        iv = ix_['X1,6']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X2,6']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C16',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C16')
        iv = ix_['X3,5']
        pbm.A[ig,iv] = float(-4.0)+pbm.A[ig,iv]
        iv = ix_['X1,7']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X2,7']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
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
        pbm.gconst = arrset(pbm.gconst,ig_['C1'],float(0.652981))
        pbm.gconst = arrset(pbm.gconst,ig_['C2'],float(0.281941))
        pbm.gconst = arrset(pbm.gconst,ig_['C3'],float(3.705233))
        pbm.gconst = arrset(pbm.gconst,ig_['C4'],float(47.00022))
        pbm.gconst = arrset(pbm.gconst,ig_['C5'],float(47.02972))
        pbm.gconst = arrset(pbm.gconst,ig_['C6'],float(0.08005))
        pbm.gconst = arrset(pbm.gconst,ig_['C7'],float(0.08813))
        pbm.gconst = arrset(pbm.gconst,ig_['C8'],float(0.04829))
        pbm.gconst = arrset(pbm.gconst,ig_['C11'],float(0.0022725))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),1.e-12)
        pb.xupper = np.full((pb.n,1),+float('inf'))
        pb.xlower[ix_['X13,2']] = 0.0155
        pb.xupper[ix_['X13,2']] = 0.0155
        pb.xlower[ix_['X13,3']] = 0.0211275
        pb.xupper[ix_['X13,3']] = 0.0211275
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(0.1))
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
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for K in range(int(v_['1']),int(v_['NSETS'])+1):
            v_['RNSK'] = v_['NS'+str(K)]
            v_['NSK'] = int(np.fix(v_['RNSK']))
            for J in range(int(v_['1']),int(v_['NSK'])+1):
                ename = 'A'+str(J)+','+str(K)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eXLOGX')
                ielftype = arrset(ielftype, ie, iet_["eXLOGX"])
                vname = 'X'+str(J)+','+str(K)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
                posev = find(elftv[ielftype[ie]],lambda x:x=='X')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        v_['K'] = 1
        for J in range(int(v_['1']),int(v_['4'])+1):
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eXLOGX4')
            ielftype = arrset(ielftype, ie, iet_["eXLOGX4"])
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(J)+','+str(int(v_['K']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X1,1'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X2,1'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X3,1'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X4,1'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y4')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        v_['K'] = 2
        for J in range(int(v_['1']),int(v_['13'])+1):
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eXLOGX13')
            ielftype = arrset(ielftype, ie, iet_["eXLOGX13"])
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(J)+','+str(int(v_['K']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X1,2'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X2,2'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X3,2'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X4,2'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y4')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X5,2'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y5')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X6,2'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y6')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X7,2'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y7')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X8,2'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y8')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X9,2'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y9')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X10,2'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y10')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X11,2'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y11')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X12,2'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y12')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X13,2'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y13')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        v_['K'] = 3
        for J in range(int(v_['1']),int(v_['18'])+1):
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eXLOGX18')
            ielftype = arrset(ielftype, ie, iet_["eXLOGX18"])
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(J)+','+str(int(v_['K']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X1,3'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X2,3'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X3,3'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X4,3'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y4')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X5,3'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y5')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X6,3'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y6')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X7,3'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y7')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X8,3'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y8')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X9,3'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y9')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X10,3'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y10')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X11,3'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y11')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X12,3'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y12')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X13,3'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y13')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X14,3'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y14')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X15,3'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y15')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X16,3'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y16')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X17,3'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y17')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X18,3'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y18')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        v_['K'] = 4
        for J in range(int(v_['1']),int(v_['3'])+1):
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eXLOGX3')
            ielftype = arrset(ielftype, ie, iet_["eXLOGX3"])
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(J)+','+str(int(v_['K']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X1,4'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X2,4'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X3,4'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        v_['K'] = 5
        for J in range(int(v_['1']),int(v_['3'])+1):
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eXLOGX3')
            ielftype = arrset(ielftype, ie, iet_["eXLOGX3"])
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(J)+','+str(int(v_['K']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X1,5'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X2,5'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X3,5'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        v_['K'] = 6
        for J in range(int(v_['1']),int(v_['2'])+1):
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eXLOGX2')
            ielftype = arrset(ielftype, ie, iet_["eXLOGX2"])
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(J)+','+str(int(v_['K']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X1,6'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X2,6'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        v_['K'] = 7
        for J in range(int(v_['1']),int(v_['2'])+1):
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eXLOGX2')
            ielftype = arrset(ielftype, ie, iet_["eXLOGX2"])
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(J)+','+str(int(v_['K']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X1,7'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(J)+','+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X2,7'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.e-12,None,0.1)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for K in range(int(v_['1']),int(v_['NSETS'])+1):
            v_['RNSK'] = v_['NS'+str(K)]
            v_['NSK'] = int(np.fix(v_['RNSK']))
            for J in range(int(v_['1']),int(v_['NSK'])+1):
                ig = ig_['OBJ']
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['A'+str(J)+','+str(K)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
                posel = posel+1
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['B'+str(J)+','+str(K)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               -1910.344724
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
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
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "OLR2-MY-45-14"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eXLOGX(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        LOGX = np.log(EV_[0])
        f_   = EV_[0]*LOGX
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = LOGX+1.0
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 1.0/EV_[0]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eXLOGX2(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((2,3))
        IV_ = np.zeros(2)
        U_[1,0] = U_[1,0]+1
        U_[1,1] = U_[1,1]+1
        U_[0,2] = U_[0,2]+1
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        LOGX = np.log(IV_[1])
        f_   = IV_[0]*LOGX
        if not isinstance( f_, float ):
            f_   = f_.item();
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
    def eXLOGX3(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((2,4))
        IV_ = np.zeros(2)
        U_[1,0] = U_[1,0]+1
        U_[1,1] = U_[1,1]+1
        U_[1,2] = U_[1,2]+1
        U_[0,3] = U_[0,3]+1
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        LOGX = np.log(IV_[1])
        f_   = IV_[0]*LOGX
        if not isinstance( f_, float ):
            f_   = f_.item();
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
    def eXLOGX4(pbm,nargout,*args):

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
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        LOGX = np.log(IV_[1])
        f_   = IV_[0]*LOGX
        if not isinstance( f_, float ):
            f_   = f_.item();
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
    def eXLOGX13(pbm,nargout,*args):

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
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        LOGX = np.log(IV_[1])
        f_   = IV_[0]*LOGX
        if not isinstance( f_, float ):
            f_   = f_.item();
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
    def eXLOGX18(pbm,nargout,*args):

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
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        LOGX = np.log(IV_[1])
        f_   = IV_[0]*LOGX
        if not isinstance( f_, float ):
            f_   = f_.item();
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

