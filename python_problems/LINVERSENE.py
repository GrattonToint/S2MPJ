from s2mpjlib import *
class  LINVERSENE(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
#    Problem : LINVERSENE
#    *********
# 
#    The problem is to find the positive definite lower bidiagonal
#    matrix L such that the matrix L(inv)L(inv-transp) best approximates,
#    in the Frobenius norm, a given symmetric target matrix T.
#    More precisely, one is  interested in the positive definite lower
#    bidiagonal L such that
# 
#         || L T L(transp) - I ||     is minimum.
#                                F
# 
#    The positive definite character of L is imposed by requiring
#    that all its diagonal entries to be at least equal to EPSILON,
#    a strictly positive real number.
# 
#    Many variants of the problem can be obtained by varying the target
#    matrix T and the scalar EPSILON.  In the present problem,
#    a) T is chosen to be pentadiagonal with T(i,j) = sin(i)cos(j) (j .leq. i)
#    b) EPSILON = 1.D-8
# 
#    Source:
#    Ph. Toint, private communication, 1991.
# 
#    SIF input: Ph. Toint, March 1991.
#    Bound-constrained nonlinear equations version: Nick Gould, June 2019.
# 
#    classification = "NOR2-AN-V-V"
# 
#    Dimension of the matrix
# 
#           Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER  n = 19    original value
# IE N                   100            $-PARAMETER  n = 199
# IE N                   500            $-PARAMETER  n = 999
# IE N                   1000           $-PARAMETER  n = 1999
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'LINVERSENE'

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
        if nargin<1:
            v_['N'] = int(10);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
# IE N                   4000           $-PARAMETER  n = 1999
# IE N                   10000          $-PARAMETER  n = 19999
        v_['EPSILON'] = 1.0e-8
        v_['ROOTP5'] = np.sqrt(0.5)
        v_['1'] = 1
        v_['2'] = 2
        v_['3'] = 3
        v_['4'] = 4
        v_['N-1'] = -1+v_['N']
        v_['N-2'] = -2+v_['N']
        for J in range(int(v_['1']),int(v_['N-2'])+1):
            v_['J+2'] = 2+J
            v_['RJ'] = float(J)
            v_['COSJ'] = np.cos(v_['RJ'])
            for I in range(int(J),int(v_['J+2'])+1):
                v_['RI'] = float(I)
                v_['SINI'] = np.sin(v_['RI'])
                v_['T'+str(I)+','+str(J)] = v_['SINI']*v_['COSJ']
        v_['RN-1'] = float(v_['N-1'])
        v_['SINI'] = np.sin(v_['RN-1'])
        v_['COSJ'] = np.cos(v_['RN-1'])
        v_['T'+str(int(v_['N-1']))+','+str(int(v_['N-1']))] = v_['SINI']*v_['COSJ']
        v_['RN'] = float(v_['N'])
        v_['SINI'] = np.sin(v_['RN'])
        v_['T'+str(int(v_['N']))+','+str(int(v_['N-1']))] = v_['SINI']*v_['COSJ']
        v_['COSJ'] = np.cos(v_['RN'])
        v_['T'+str(int(v_['N']))+','+str(int(v_['N']))] = v_['SINI']*v_['COSJ']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        pb.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            [iv,ix_,_] = s2mpj_ii('A'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'A'+str(I))
            [iv,ix_,_] = s2mpj_ii('B'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'B'+str(I))
        [iv,ix_,_] = s2mpj_ii('A'+str(int(v_['N'])),ix_)
        pb.xnames=arrset(pb.xnames,iv,'A'+str(int(v_['N'])))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for J in range(int(v_['1']),int(v_['N-2'])+1):
            v_['J+1'] = 1+J
            v_['J+2'] = 2+J
            [ig,ig_,_] = s2mpj_ii('O'+str(J)+','+str(J),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'O'+str(J)+','+str(J))
            [ig,ig_,_] = s2mpj_ii('O'+str(int(v_['J+1']))+','+str(J),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'O'+str(int(v_['J+1']))+','+str(J))
            pbm.gscale = arrset(pbm.gscale,ig,float(v_['ROOTP5']))
            [ig,ig_,_] = s2mpj_ii('O'+str(int(v_['J+2']))+','+str(J),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'O'+str(int(v_['J+2']))+','+str(J))
            pbm.gscale = arrset(pbm.gscale,ig,float(v_['ROOTP5']))
        [ig,ig_,_] = s2mpj_ii('O'+str(int(v_['N-1']))+','+str(int(v_['N-1'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'O'+str(int(v_['N-1']))+','+str(int(v_['N-1'])))
        [ig,ig_,_] = s2mpj_ii('O'+str(int(v_['N']))+','+str(int(v_['N-1'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'O'+str(int(v_['N']))+','+str(int(v_['N-1'])))
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['ROOTP5']))
        [ig,ig_,_] = s2mpj_ii('O'+str(int(v_['N']))+','+str(int(v_['N'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'O'+str(int(v_['N']))+','+str(int(v_['N'])))
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
        for I in range(int(v_['1']),int(v_['N'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['O'+str(I)+','+str(I)],float(1.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower = np.zeros((pb.n,1))
        for I in range(int(v_['1']),int(v_['N'])+1):
            pb.xlower[ix_['A'+str(I)]] = v_['EPSILON']
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(-1.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        ename = 'S'+str(int(v_['1']))+','+str(int(v_['1']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['1']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'S'+str(int(v_['1']))+','+str(int(v_['1']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['1']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'S'+str(int(v_['2']))+','+str(int(v_['1']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['2']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'S'+str(int(v_['2']))+','+str(int(v_['1']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['1']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'V'+str(int(v_['2']))+','+str(int(v_['1']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'B'+str(int(v_['1']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'V'+str(int(v_['2']))+','+str(int(v_['1']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['1']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'S'+str(int(v_['3']))+','+str(int(v_['1']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['3']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'S'+str(int(v_['3']))+','+str(int(v_['1']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['1']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'V'+str(int(v_['3']))+','+str(int(v_['1']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'B'+str(int(v_['2']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'V'+str(int(v_['3']))+','+str(int(v_['1']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['1']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'S'+str(int(v_['2']))+','+str(int(v_['2']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['2']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'S'+str(int(v_['2']))+','+str(int(v_['2']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['2']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'U'+str(int(v_['2']))+','+str(int(v_['2']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['2']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'U'+str(int(v_['2']))+','+str(int(v_['2']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'B'+str(int(v_['1']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'V'+str(int(v_['2']))+','+str(int(v_['2']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['2']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'V'+str(int(v_['2']))+','+str(int(v_['2']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'B'+str(int(v_['1']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'W'+str(int(v_['2']))+','+str(int(v_['2']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'B'+str(int(v_['1']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'W'+str(int(v_['2']))+','+str(int(v_['2']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'B'+str(int(v_['1']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'S'+str(int(v_['3']))+','+str(int(v_['2']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['3']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'S'+str(int(v_['3']))+','+str(int(v_['2']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['2']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'U'+str(int(v_['3']))+','+str(int(v_['2']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['3']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'U'+str(int(v_['3']))+','+str(int(v_['2']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'B'+str(int(v_['1']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'V'+str(int(v_['3']))+','+str(int(v_['2']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['2']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'V'+str(int(v_['3']))+','+str(int(v_['2']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'B'+str(int(v_['2']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'W'+str(int(v_['3']))+','+str(int(v_['2']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'B'+str(int(v_['2']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'W'+str(int(v_['3']))+','+str(int(v_['2']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'B'+str(int(v_['1']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'S'+str(int(v_['3']))+','+str(int(v_['3']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['3']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'S'+str(int(v_['3']))+','+str(int(v_['3']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['3']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'U'+str(int(v_['3']))+','+str(int(v_['3']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['3']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'U'+str(int(v_['3']))+','+str(int(v_['3']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'B'+str(int(v_['2']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'V'+str(int(v_['3']))+','+str(int(v_['3']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['3']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'V'+str(int(v_['3']))+','+str(int(v_['3']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'B'+str(int(v_['2']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'W'+str(int(v_['3']))+','+str(int(v_['3']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'B'+str(int(v_['2']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'W'+str(int(v_['3']))+','+str(int(v_['3']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'B'+str(int(v_['2']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['4']),int(v_['N'])+1):
            v_['I-1'] = -1+I
            v_['I-2'] = -2+I
            ename = 'S'+str(I)+','+str(int(v_['I-2']))
            [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
            if newelt:
                pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
                ielftype = arrset( ielftype,ie,iet_['en2PR'])
            vname = 'A'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'S'+str(I)+','+str(int(v_['I-2']))
            [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
            if newelt:
                pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
                ielftype = arrset( ielftype,ie,iet_['en2PR'])
            vname = 'A'+str(int(v_['I-2']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'V'+str(I)+','+str(int(v_['I-2']))
            [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
            if newelt:
                pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
                ielftype = arrset( ielftype,ie,iet_['en2PR'])
            vname = 'B'+str(int(v_['I-1']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'V'+str(I)+','+str(int(v_['I-2']))
            [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
            if newelt:
                pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
                ielftype = arrset( ielftype,ie,iet_['en2PR'])
            vname = 'A'+str(int(v_['I-2']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            for J in range(int(v_['I-1']),int(I)+1):
                v_['J-1'] = -1+J
                ename = 'S'+str(I)+','+str(J)
                [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
                if newelt:
                    pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
                    ielftype = arrset( ielftype,ie,iet_['en2PR'])
                vname = 'A'+str(I)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
                posev = find(elftv[ielftype[ie]],lambda x:x=='X')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'A'+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
                posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                ename = 'U'+str(I)+','+str(J)
                [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
                if newelt:
                    pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
                    ielftype = arrset( ielftype,ie,iet_['en2PR'])
                vname = 'A'+str(I)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
                posev = find(elftv[ielftype[ie]],lambda x:x=='X')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'B'+str(int(v_['J-1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
                posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                ename = 'V'+str(I)+','+str(J)
                [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
                if newelt:
                    pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
                    ielftype = arrset( ielftype,ie,iet_['en2PR'])
                vname = 'B'+str(int(v_['I-1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
                posev = find(elftv[ielftype[ie]],lambda x:x=='X')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'A'+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
                posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                ename = 'W'+str(I)+','+str(J)
                [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
                if newelt:
                    pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
                    ielftype = arrset( ielftype,ie,iet_['en2PR'])
                vname = 'B'+str(int(v_['I-1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
                posev = find(elftv[ielftype[ie]],lambda x:x=='X')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'B'+str(int(v_['J-1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,-1.0)
                posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        ig = ig_['O'+str(int(v_['1']))+','+str(int(v_['1']))]
        posel = len(pbm.grelt[ig])
        pbm.grelt  = (
              loaset(pbm.grelt,ig,posel,ie_['S'+str(int(v_['1']))+','+str(int(v_['1']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw  = (
              loaset(pbm.grelw,ig,posel,float(v_['T'+str(int(v_['1']))+','+str(int(v_['1']))])))
        ig = ig_['O'+str(int(v_['2']))+','+str(int(v_['1']))]
        posel = len(pbm.grelt[ig])
        pbm.grelt  = (
              loaset(pbm.grelt,ig,posel,ie_['S'+str(int(v_['2']))+','+str(int(v_['1']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw  = (
              loaset(pbm.grelw,ig,posel,float(v_['T'+str(int(v_['2']))+','+str(int(v_['1']))])))
        posel = len(pbm.grelt[ig])
        pbm.grelt  = (
              loaset(pbm.grelt,ig,posel,ie_['V'+str(int(v_['2']))+','+str(int(v_['1']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw  = (
              loaset(pbm.grelw,ig,posel,float(v_['T'+str(int(v_['1']))+','+str(int(v_['1']))])))
        ig = ig_['O'+str(int(v_['3']))+','+str(int(v_['1']))]
        posel = len(pbm.grelt[ig])
        pbm.grelt  = (
              loaset(pbm.grelt,ig,posel,ie_['S'+str(int(v_['3']))+','+str(int(v_['1']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw  = (
              loaset(pbm.grelw,ig,posel,float(v_['T'+str(int(v_['3']))+','+str(int(v_['1']))])))
        posel = len(pbm.grelt[ig])
        pbm.grelt  = (
              loaset(pbm.grelt,ig,posel,ie_['V'+str(int(v_['3']))+','+str(int(v_['1']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw  = (
              loaset(pbm.grelw,ig,posel,float(v_['T'+str(int(v_['2']))+','+str(int(v_['1']))])))
        ig = ig_['O'+str(int(v_['2']))+','+str(int(v_['2']))]
        posel = len(pbm.grelt[ig])
        pbm.grelt  = (
              loaset(pbm.grelt,ig,posel,ie_['S'+str(int(v_['2']))+','+str(int(v_['2']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw  = (
              loaset(pbm.grelw,ig,posel,float(v_['T'+str(int(v_['2']))+','+str(int(v_['2']))])))
        posel = len(pbm.grelt[ig])
        pbm.grelt  = (
              loaset(pbm.grelt,ig,posel,ie_['U'+str(int(v_['2']))+','+str(int(v_['2']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw  = (
              loaset(pbm.grelw,ig,posel,float(v_['T'+str(int(v_['2']))+','+str(int(v_['1']))])))
        posel = len(pbm.grelt[ig])
        pbm.grelt  = (
              loaset(pbm.grelt,ig,posel,ie_['V'+str(int(v_['2']))+','+str(int(v_['2']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw  = (
              loaset(pbm.grelw,ig,posel,float(v_['T'+str(int(v_['2']))+','+str(int(v_['1']))])))
        posel = len(pbm.grelt[ig])
        pbm.grelt  = (
              loaset(pbm.grelt,ig,posel,ie_['W'+str(int(v_['2']))+','+str(int(v_['2']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw  = (
              loaset(pbm.grelw,ig,posel,float(v_['T'+str(int(v_['1']))+','+str(int(v_['1']))])))
        ig = ig_['O'+str(int(v_['3']))+','+str(int(v_['2']))]
        posel = len(pbm.grelt[ig])
        pbm.grelt  = (
              loaset(pbm.grelt,ig,posel,ie_['S'+str(int(v_['3']))+','+str(int(v_['2']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw  = (
              loaset(pbm.grelw,ig,posel,float(v_['T'+str(int(v_['3']))+','+str(int(v_['2']))])))
        posel = len(pbm.grelt[ig])
        pbm.grelt  = (
              loaset(pbm.grelt,ig,posel,ie_['U'+str(int(v_['3']))+','+str(int(v_['2']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw  = (
              loaset(pbm.grelw,ig,posel,float(v_['T'+str(int(v_['3']))+','+str(int(v_['1']))])))
        posel = len(pbm.grelt[ig])
        pbm.grelt  = (
              loaset(pbm.grelt,ig,posel,ie_['V'+str(int(v_['3']))+','+str(int(v_['2']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw  = (
              loaset(pbm.grelw,ig,posel,float(v_['T'+str(int(v_['2']))+','+str(int(v_['2']))])))
        posel = len(pbm.grelt[ig])
        pbm.grelt  = (
              loaset(pbm.grelt,ig,posel,ie_['W'+str(int(v_['3']))+','+str(int(v_['2']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw  = (
              loaset(pbm.grelw,ig,posel,float(v_['T'+str(int(v_['2']))+','+str(int(v_['1']))])))
        ig = ig_['O'+str(int(v_['3']))+','+str(int(v_['3']))]
        posel = len(pbm.grelt[ig])
        pbm.grelt  = (
              loaset(pbm.grelt,ig,posel,ie_['S'+str(int(v_['3']))+','+str(int(v_['3']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw  = (
              loaset(pbm.grelw,ig,posel,float(v_['T'+str(int(v_['3']))+','+str(int(v_['3']))])))
        posel = len(pbm.grelt[ig])
        pbm.grelt  = (
              loaset(pbm.grelt,ig,posel,ie_['U'+str(int(v_['3']))+','+str(int(v_['3']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw  = (
              loaset(pbm.grelw,ig,posel,float(v_['T'+str(int(v_['3']))+','+str(int(v_['2']))])))
        posel = len(pbm.grelt[ig])
        pbm.grelt  = (
              loaset(pbm.grelt,ig,posel,ie_['V'+str(int(v_['3']))+','+str(int(v_['3']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw  = (
              loaset(pbm.grelw,ig,posel,float(v_['T'+str(int(v_['3']))+','+str(int(v_['2']))])))
        posel = len(pbm.grelt[ig])
        pbm.grelt  = (
              loaset(pbm.grelt,ig,posel,ie_['W'+str(int(v_['3']))+','+str(int(v_['3']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw  = (
              loaset(pbm.grelw,ig,posel,float(v_['T'+str(int(v_['2']))+','+str(int(v_['2']))])))
        for I in range(int(v_['4']),int(v_['N'])+1):
            v_['I-1'] = -1+I
            v_['I-2'] = -2+I
            ig = ig_['O'+str(I)+','+str(int(v_['I-2']))]
            posel = len(pbm.grelt[ig])
            pbm.grelt  = (
                  loaset(pbm.grelt,ig,posel,ie_['S'+str(I)+','+str(int(v_['I-2']))]))
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw  = (
                  loaset(pbm.grelw,ig,posel,float(v_['T'+str(I)+','+str(int(v_['I-2']))])))
            posel = len(pbm.grelt[ig])
            pbm.grelt  = (
                  loaset(pbm.grelt,ig,posel,ie_['V'+str(I)+','+str(int(v_['I-2']))]))
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw  = (
                  loaset(pbm.grelw,ig,posel,float(v_['T'+str(int(v_['I-1']))+','+str(int(v_['I-2']))])))
            ig = ig_['O'+str(I)+','+str(int(v_['I-1']))]
            posel = len(pbm.grelt[ig])
            pbm.grelt  = (
                  loaset(pbm.grelt,ig,posel,ie_['S'+str(I)+','+str(int(v_['I-1']))]))
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw  = (
                  loaset(pbm.grelw,ig,posel,float(v_['T'+str(I)+','+str(int(v_['I-1']))])))
            posel = len(pbm.grelt[ig])
            pbm.grelt  = (
                  loaset(pbm.grelt,ig,posel,ie_['U'+str(I)+','+str(int(v_['I-1']))]))
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw  = (
                  loaset(pbm.grelw,ig,posel,float(v_['T'+str(I)+','+str(int(v_['I-2']))])))
            posel = len(pbm.grelt[ig])
            pbm.grelt  = (
                  loaset(pbm.grelt,ig,posel,ie_['V'+str(I)+','+str(int(v_['I-1']))]))
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw  = (
                  loaset(pbm.grelw,ig,posel,float(v_['T'+str(int(v_['I-1']))+','+str(int(v_['I-1']))])))
            posel = len(pbm.grelt[ig])
            pbm.grelt  = (
                  loaset(pbm.grelt,ig,posel,ie_['W'+str(I)+','+str(int(v_['I-1']))]))
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw  = (
                  loaset(pbm.grelw,ig,posel,float(v_['T'+str(int(v_['I-1']))+','+str(int(v_['I-2']))])))
            ig = ig_['O'+str(I)+','+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['S'+str(I)+','+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['T'+str(I)+','+str(I)]))
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['U'+str(I)+','+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw  = (
                  loaset(pbm.grelw,ig,posel,float(v_['T'+str(I)+','+str(int(v_['I-1']))])))
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['V'+str(I)+','+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw  = (
                  loaset(pbm.grelw,ig,posel,float(v_['T'+str(I)+','+str(int(v_['I-1']))])))
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['W'+str(I)+','+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw  = (
                  loaset(pbm.grelw,ig,posel,float(v_['T'+str(int(v_['I-1']))+','+str(int(v_['I-1']))])))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN(10)           6.00000000
# LO SOLTN(100)          68.0000000
# LO SOLTN(500)          340.000000
# LO SOLTN(1000)         ???
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "NOR2-AN-V-V"
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

