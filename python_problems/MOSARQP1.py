from s2xlib import *
class  MOSARQP1(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : MOSARQP1
#    *********
# 
#    A convex quadratic problem, with variable dimensions.
#    In this problem, half the linear constraints are active at the solution.
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
#    classification = "QLR2-AN-V-V"
# 
#    Problem variants: these are distinguished by the triplet ( N, M, COND ),
#    where: - N (nb of variables) must be even and have an integer square root
#           - M (nb of constraints) must be at least sqrt(N) 
#             and at most N - sqrt(N)
#           - COND (problem conditioning) is a positive integer
#    Except for the first, the instances suggested are those used by Morales
#    and Sargent.
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'MOSARQP1'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'MOSARQP1'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(36);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
        if nargin<2:
            v_['M'] = int(10);  #  SIF file default value
        else:
            v_['M'] = int(args[1])
        if nargin<3:
            v_['COND'] = float(2.0);  #  SIF file default value
        else:
            v_['COND'] = float(args[2])
#           Alternative values for the SIF file parameters:
# IE N                   100            $-PARAMETER     original value
# IE M                   10             $-PARAMETER     original value
# RE COND                3.0            $-PARAMETER     original value
# IE N                   900            $-PARAMETER
# IE M                   30             $-PARAMETER
# RE COND                1.0            $-PARAMETER
# IE N                   900            $-PARAMETER
# IE M                   30             $-PARAMETER
# RE COND                2.0            $-PARAMETER
# IE N                   900            $-PARAMETER
# IE M                   30             $-PARAMETER
# RE COND                3.0            $-PARAMETER
# IE N                   900            $-PARAMETER
# IE M                   60             $-PARAMETER
# RE COND                1.0            $-PARAMETER
# IE N                   900            $-PARAMETER
# IE M                   60             $-PARAMETER
# RE COND                2.0            $-PARAMETER
# IE N                   900            $-PARAMETER
# IE M                   60             $-PARAMETER
# RE COND                3.0            $-PARAMETER
# IE N                   900            $-PARAMETER
# IE M                   90             $-PARAMETER
# RE COND                1.0            $-PARAMETER
# IE N                   900            $-PARAMETER
# IE M                   90             $-PARAMETER
# RE COND                2.0            $-PARAMETER
# IE N                   900            $-PARAMETER
# IE M                   90             $-PARAMETER
# RE COND                3.0            $-PARAMETER
# IE N                   900            $-PARAMETER
# IE M                   120            $-PARAMETER
# RE COND                1.0            $-PARAMETER
# IE N                   900            $-PARAMETER
# IE M                   120            $-PARAMETER
# RE COND                2.0            $-PARAMETER
# IE N                   900            $-PARAMETER
# IE M                   120            $-PARAMETER
# RE COND                3.0            $-PARAMETER
# IE N                   900            $-PARAMETER
# IE M                   300            $-PARAMETER
# RE COND                1.0            $-PARAMETER
# IE N                   900            $-PARAMETER
# IE M                   300            $-PARAMETER
# RE COND                2.0            $-PARAMETER
# IE N                   900            $-PARAMETER
# IE M                   300            $-PARAMETER
# RE COND                3.0            $-PARAMETER
# IE N                   900            $-PARAMETER
# IE M                   600            $-PARAMETER
# RE COND                1.0            $-PARAMETER
# IE N                   900            $-PARAMETER
# IE M                   600            $-PARAMETER
# RE COND                2.0            $-PARAMETER
# IE N                   900            $-PARAMETER
# IE M                   600            $-PARAMETER
# RE COND                3.0            $-PARAMETER
# IE N                   2500           $-PARAMETER
# IE M                   700            $-PARAMETER
# RE COND                1.0            $-PARAMETER
# IE N                   2500           $-PARAMETER
# IE M                   700            $-PARAMETER
# RE COND                2.0            $-PARAMETER
# IE N                   2500           $-PARAMETER
# IE M                   700            $-PARAMETER
# RE COND                3.0            $-PARAMETER
        v_['1'] = 1
        v_['2'] = 2
        v_['N-1'] = -1+v_['N']
        v_['RN-1'] = float(v_['N-1'])
        v_['RN'] = float(v_['N'])
        v_['M-1'] = -1+v_['M']
        v_['RNP'] = 0.1+v_['RN']
        v_['RRTN'] = np.sqrt(v_['RNP'])
        v_['RTN'] = int(np.fix(v_['RRTN']))
        v_['RTN-1'] = -1+v_['RTN']
        v_['RTN-2'] = -2+v_['RTN']
        v_['RTN+1'] = 1+v_['RTN']
        v_['2RTN'] = v_['RTN']+v_['RTN']
        v_['M-RTN+1'] = v_['M']-v_['RTN-1']
        for I in range(int(v_['1']),int(v_['N-1'])+1,int(v_['2'])):
            v_['I+1'] = 1+I
            v_['XC'+str(I)] = -1.0
            v_['XC'+str(int(v_['I+1']))] = 1.0
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
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2x_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [ig,ig_,_] = s2x_ii('OBJ',ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = float(v_['C'+str(I)])+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('CS'+str(int(v_['1'])),ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'CS'+str(int(v_['1'])))
        iv = ix_['X'+str(int(v_['1']))]
        pbm.A[ig,iv] = float(4.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('CS'+str(int(v_['1'])),ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'CS'+str(int(v_['1'])))
        iv = ix_['X'+str(int(v_['RTN+1']))]
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X'+str(int(v_['2']))]
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        for I in range(int(v_['2']),int(v_['RTN-1'])+1):
            v_['I+1'] = 1+I
            v_['I-1'] = -1+I
            v_['I+RTN'] = I+v_['RTN']
            [ig,ig_,_] = s2x_ii('CS'+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'CS'+str(I))
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = float(4.0)+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['I+RTN']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['I-1']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['I+1']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('CS'+str(int(v_['RTN'])),ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'CS'+str(int(v_['RTN'])))
        iv = ix_['X'+str(int(v_['RTN']))]
        pbm.A[ig,iv] = float(4.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('CS'+str(int(v_['RTN'])),ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'CS'+str(int(v_['RTN'])))
        iv = ix_['X'+str(int(v_['RTN-1']))]
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X'+str(int(v_['2RTN']))]
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        v_['JS'] = v_['RTN']
        for J in range(int(v_['RTN+1']),int(v_['M-RTN+1'])+1,int(v_['RTN'])):
            v_['J+1'] = 1+J
            v_['JS'] = J+v_['RTN-1']
            v_['JS-1'] = -1+v_['JS']
            v_['J-RTN'] = J-v_['RTN']
            v_['J+RTN'] = J+v_['RTN']
            [ig,ig_,_] = s2x_ii('CS'+str(J),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'CS'+str(J))
            iv = ix_['X'+str(J)]
            pbm.A[ig,iv] = float(4.0)+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['J+1']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['J-RTN']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['J+RTN']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            for I in range(int(v_['J+1']),int(v_['JS-1'])+1):
                v_['I+1'] = 1+I
                v_['I-1'] = -1+I
                v_['I+RTN'] = I+v_['RTN']
                v_['I-RTN'] = I-v_['RTN']
                [ig,ig_,_] = s2x_ii('CS'+str(I),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'CS'+str(I))
                iv = ix_['X'+str(I)]
                pbm.A[ig,iv] = float(4.0)+pbm.A[ig,iv]
                iv = ix_['X'+str(int(v_['I-1']))]
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
                iv = ix_['X'+str(int(v_['I+1']))]
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
                iv = ix_['X'+str(int(v_['I-RTN']))]
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
                iv = ix_['X'+str(int(v_['I+RTN']))]
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            v_['JS+RTN'] = v_['JS']+v_['RTN']
            v_['JS-RTN'] = v_['JS']-v_['RTN']
            [ig,ig_,_] = s2x_ii('CS'+str(int(v_['JS'])),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'CS'+str(int(v_['JS'])))
            iv = ix_['X'+str(int(v_['JS']))]
            pbm.A[ig,iv] = float(4.0)+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['JS-1']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('CS'+str(int(v_['JS'])),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'CS'+str(int(v_['JS'])))
            iv = ix_['X'+str(int(v_['JS-RTN']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['JS+RTN']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        v_['K'] = 1+v_['JS']
        for I in range(int(v_['K']),int(v_['M'])+1,int(v_['M'])):
            v_['K+1'] = 1+v_['K']
            v_['K+RTN'] = v_['K']+v_['RTN']
            v_['K-RTN'] = v_['K']-v_['RTN']
            [ig,ig_,_] = s2x_ii('CS'+str(int(v_['K'])),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'CS'+str(int(v_['K'])))
            iv = ix_['X'+str(int(v_['K']))]
            pbm.A[ig,iv] = float(4.0)+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['K+1']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('CS'+str(int(v_['K'])),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'CS'+str(int(v_['K'])))
            iv = ix_['X'+str(int(v_['K-RTN']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['K+RTN']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        v_['K'] = 1+v_['K']
        for I in range(int(v_['K']),int(v_['M'])+1):
            v_['I+1'] = 1+I
            v_['I-1'] = -1+I
            v_['I+RTN'] = I+v_['RTN']
            v_['I-RTN'] = I-v_['RTN']
            [ig,ig_,_] = s2x_ii('CS'+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'CS'+str(I))
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = float(4.0)+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['I-1']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['I+1']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['I-RTN']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['I+RTN']))]
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
        pbm.gconst = arrset(pbm.gconst,ig_['CS'+str(int(v_['1']))],float(0.5))
        pbm.gconst = arrset(pbm.gconst,ig_['CS'+str(int(v_['RTN']))],float(0.5))
        v_['K'] = v_['RTN+1']
        for J in range(int(v_['RTN+1']),int(v_['M-RTN+1'])+1,int(v_['RTN'])):
            v_['K'] = 1+v_['K']
            for I in range(int(v_['1']),int(v_['RTN-2'])+1):
                v_['K'] = 1+v_['K']
                pbm.gconst = arrset(pbm.gconst,ig_['CS'+str(int(v_['K']))],float(-0.5))
            v_['K'] = 1+v_['K']
        v_['K'] = 1+v_['K']
        for J in range(int(v_['K']),int(v_['M'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['CS'+str(J)],float(-0.5))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(0.5))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'X')
        [it,iet_,_] = s2x_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['1']),int(v_['N'])+1):
            ename = 'XSQ'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
            ielftype = arrset(ielftype, ie, iet_["eSQ"])
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,0.5)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['1']),int(v_['NNZ'])+1):
            v_['RKI'] = v_['K'+str(I)]
            v_['KI'] = int(np.fix(v_['RKI']))
            v_['I-1'] = -1+I
            for J in range(int(v_['1']),int(v_['I-1'])+1):
                v_['RKJ'] = v_['K'+str(J)]
                v_['KJ'] = int(np.fix(v_['RKJ']))
                ename = 'P'+str(I)+','+str(J)
                [ie,ie_,_] = s2x_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
                ielftype = arrset(ielftype, ie, iet_["en2PR"])
                vname = 'X'+str(int(v_['KI']))
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,0.5)
                posev = find(elftv[ielftype[ie]],lambda x:x=='X')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(int(v_['KJ']))
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,0.5)
                posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['TMP'] = 0.5*v_['D'+str(I)]
            ig = ig_['OBJ']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['XSQ'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['TMP']))
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
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['WIJ']))
            v_['TMP'] = v_['DY'+str(I)]*v_['Y'+str(I)]
            v_['WII'] = v_['TMP']*v_['-2/YN2']
            v_['TMP'] = v_['Y'+str(I)]*v_['Y'+str(I)]
            v_['TMP'] = v_['TMP']*v_['DD/2']
            v_['WII'] = v_['WII']+v_['TMP']
            ig = ig_['OBJ']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['XSQ'+str(int(v_['KI']))])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['WII']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle+pb.neq,pb.m)] = np.zeros((pb.nge,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "QLR2-AN-V-V"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSQ(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[0]+EV_[0]
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
