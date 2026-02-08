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
#    classification = "C-CNOR2-AN-V-V"
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
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'LINVERSENE'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
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
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            [iv,ix_,_] = s2mpj_ii('A'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'A'+str(I))
            [iv,ix_,_] = s2mpj_ii('B'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'B'+str(I))
        [iv,ix_,_] = s2mpj_ii('A'+str(int(v_['N'])),ix_)
        self.xnames=arrset(self.xnames,iv,'A'+str(int(v_['N'])))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for J in range(int(v_['1']),int(v_['N-2'])+1):
            v_['J+1'] = 1+J
            v_['J+2'] = 2+J
            [ig,ig_,_] = s2mpj_ii('O'+str(J)+','+str(J),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'O'+str(J)+','+str(J))
            [ig,ig_,_] = s2mpj_ii('O'+str(int(v_['J+1']))+','+str(J),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'O'+str(int(v_['J+1']))+','+str(J))
            self.gscale = arrset(self.gscale,ig,float(v_['ROOTP5']))
            [ig,ig_,_] = s2mpj_ii('O'+str(int(v_['J+2']))+','+str(J),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'O'+str(int(v_['J+2']))+','+str(J))
            self.gscale = arrset(self.gscale,ig,float(v_['ROOTP5']))
        [ig,ig_,_] = s2mpj_ii('O'+str(int(v_['N-1']))+','+str(int(v_['N-1'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'O'+str(int(v_['N-1']))+','+str(int(v_['N-1'])))
        [ig,ig_,_] = s2mpj_ii('O'+str(int(v_['N']))+','+str(int(v_['N-1'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'O'+str(int(v_['N']))+','+str(int(v_['N-1'])))
        self.gscale = arrset(self.gscale,ig,float(v_['ROOTP5']))
        [ig,ig_,_] = s2mpj_ii('O'+str(int(v_['N']))+','+str(int(v_['N'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'O'+str(int(v_['N']))+','+str(int(v_['N'])))
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
            self.gconst = arrset(self.gconst,ig_['O'+str(I)+','+str(I)],float(1.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        for I in range(int(v_['1']),int(v_['N'])+1):
            self.xlower[ix_['A'+str(I)]] = v_['EPSILON']
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(-1.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        ename = 'S'+str(int(v_['1']))+','+str(int(v_['1']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['1']))
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'S'+str(int(v_['1']))+','+str(int(v_['1']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['1']))
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'S'+str(int(v_['2']))+','+str(int(v_['1']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['2']))
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'S'+str(int(v_['2']))+','+str(int(v_['1']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['1']))
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'V'+str(int(v_['2']))+','+str(int(v_['1']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'B'+str(int(v_['1']))
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'V'+str(int(v_['2']))+','+str(int(v_['1']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['1']))
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'S'+str(int(v_['3']))+','+str(int(v_['1']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['3']))
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'S'+str(int(v_['3']))+','+str(int(v_['1']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['1']))
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'V'+str(int(v_['3']))+','+str(int(v_['1']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'B'+str(int(v_['2']))
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'V'+str(int(v_['3']))+','+str(int(v_['1']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['1']))
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'S'+str(int(v_['2']))+','+str(int(v_['2']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['2']))
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'S'+str(int(v_['2']))+','+str(int(v_['2']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['2']))
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'U'+str(int(v_['2']))+','+str(int(v_['2']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['2']))
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'U'+str(int(v_['2']))+','+str(int(v_['2']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'B'+str(int(v_['1']))
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'V'+str(int(v_['2']))+','+str(int(v_['2']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['2']))
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'V'+str(int(v_['2']))+','+str(int(v_['2']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'B'+str(int(v_['1']))
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'W'+str(int(v_['2']))+','+str(int(v_['2']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'B'+str(int(v_['1']))
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'W'+str(int(v_['2']))+','+str(int(v_['2']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'B'+str(int(v_['1']))
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'S'+str(int(v_['3']))+','+str(int(v_['2']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['3']))
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'S'+str(int(v_['3']))+','+str(int(v_['2']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['2']))
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'U'+str(int(v_['3']))+','+str(int(v_['2']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['3']))
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'U'+str(int(v_['3']))+','+str(int(v_['2']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'B'+str(int(v_['1']))
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'V'+str(int(v_['3']))+','+str(int(v_['2']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['2']))
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'V'+str(int(v_['3']))+','+str(int(v_['2']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'B'+str(int(v_['2']))
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'W'+str(int(v_['3']))+','+str(int(v_['2']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'B'+str(int(v_['2']))
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'W'+str(int(v_['3']))+','+str(int(v_['2']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'B'+str(int(v_['1']))
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'S'+str(int(v_['3']))+','+str(int(v_['3']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['3']))
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'S'+str(int(v_['3']))+','+str(int(v_['3']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['3']))
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'U'+str(int(v_['3']))+','+str(int(v_['3']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['3']))
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'U'+str(int(v_['3']))+','+str(int(v_['3']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'B'+str(int(v_['2']))
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'V'+str(int(v_['3']))+','+str(int(v_['3']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'A'+str(int(v_['3']))
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'V'+str(int(v_['3']))+','+str(int(v_['3']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'B'+str(int(v_['2']))
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'W'+str(int(v_['3']))+','+str(int(v_['3']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'B'+str(int(v_['2']))
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'W'+str(int(v_['3']))+','+str(int(v_['3']))
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'B'+str(int(v_['2']))
        [iv,ix_]  = (
              s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        for I in range(int(v_['4']),int(v_['N'])+1):
            v_['I-1'] = -1+I
            v_['I-2'] = -2+I
            ename = 'S'+str(I)+','+str(int(v_['I-2']))
            [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
            if newelt:
                self.elftype = arrset(self.elftype,ie,'en2PR')
                ielftype = arrset(ielftype,ie,iet_['en2PR'])
            vname = 'A'+str(I)
            [iv,ix_]  = (
                  s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'S'+str(I)+','+str(int(v_['I-2']))
            [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
            if newelt:
                self.elftype = arrset(self.elftype,ie,'en2PR')
                ielftype = arrset(ielftype,ie,iet_['en2PR'])
            vname = 'A'+str(int(v_['I-2']))
            [iv,ix_]  = (
                  s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
            posev = np.where(elftv[ielftype[ie]]=='Y')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'V'+str(I)+','+str(int(v_['I-2']))
            [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
            if newelt:
                self.elftype = arrset(self.elftype,ie,'en2PR')
                ielftype = arrset(ielftype,ie,iet_['en2PR'])
            vname = 'B'+str(int(v_['I-1']))
            [iv,ix_]  = (
                  s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'V'+str(I)+','+str(int(v_['I-2']))
            [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
            if newelt:
                self.elftype = arrset(self.elftype,ie,'en2PR')
                ielftype = arrset(ielftype,ie,iet_['en2PR'])
            vname = 'A'+str(int(v_['I-2']))
            [iv,ix_]  = (
                  s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
            posev = np.where(elftv[ielftype[ie]]=='Y')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            for J in range(int(v_['I-1']),int(I)+1):
                v_['J-1'] = -1+J
                ename = 'S'+str(I)+','+str(J)
                [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
                if newelt:
                    self.elftype = arrset(self.elftype,ie,'en2PR')
                    ielftype = arrset(ielftype,ie,iet_['en2PR'])
                vname = 'A'+str(I)
                [iv,ix_]  = (
                      s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
                posev = np.where(elftv[ielftype[ie]]=='X')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'A'+str(J)
                [iv,ix_]  = (
                      s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
                posev = np.where(elftv[ielftype[ie]]=='Y')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                ename = 'U'+str(I)+','+str(J)
                [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
                if newelt:
                    self.elftype = arrset(self.elftype,ie,'en2PR')
                    ielftype = arrset(ielftype,ie,iet_['en2PR'])
                vname = 'A'+str(I)
                [iv,ix_]  = (
                      s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
                posev = np.where(elftv[ielftype[ie]]=='X')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'B'+str(int(v_['J-1']))
                [iv,ix_]  = (
                      s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
                posev = np.where(elftv[ielftype[ie]]=='Y')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                ename = 'V'+str(I)+','+str(J)
                [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
                if newelt:
                    self.elftype = arrset(self.elftype,ie,'en2PR')
                    ielftype = arrset(ielftype,ie,iet_['en2PR'])
                vname = 'B'+str(int(v_['I-1']))
                [iv,ix_]  = (
                      s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
                posev = np.where(elftv[ielftype[ie]]=='X')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'A'+str(J)
                [iv,ix_]  = (
                      s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
                posev = np.where(elftv[ielftype[ie]]=='Y')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                ename = 'W'+str(I)+','+str(J)
                [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
                if newelt:
                    self.elftype = arrset(self.elftype,ie,'en2PR')
                    ielftype = arrset(ielftype,ie,iet_['en2PR'])
                vname = 'B'+str(int(v_['I-1']))
                [iv,ix_]  = (
                      s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
                posev = np.where(elftv[ielftype[ie]]=='X')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'B'+str(int(v_['J-1']))
                [iv,ix_]  = (
                      s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(-1.0)))
                posev = np.where(elftv[ielftype[ie]]=='Y')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        ig = ig_['O'+str(int(v_['1']))+','+str(int(v_['1']))]
        posel = len(self.grelt[ig])
        self.grelt  = (
              loaset(self.grelt,ig,posel,ie_['S'+str(int(v_['1']))+','+str(int(v_['1']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw  = (
              loaset(self.grelw,ig,posel,float(v_['T'+str(int(v_['1']))+','+str(int(v_['1']))])))
        ig = ig_['O'+str(int(v_['2']))+','+str(int(v_['1']))]
        posel = len(self.grelt[ig])
        self.grelt  = (
              loaset(self.grelt,ig,posel,ie_['S'+str(int(v_['2']))+','+str(int(v_['1']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw  = (
              loaset(self.grelw,ig,posel,float(v_['T'+str(int(v_['2']))+','+str(int(v_['1']))])))
        posel = len(self.grelt[ig])
        self.grelt  = (
              loaset(self.grelt,ig,posel,ie_['V'+str(int(v_['2']))+','+str(int(v_['1']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw  = (
              loaset(self.grelw,ig,posel,float(v_['T'+str(int(v_['1']))+','+str(int(v_['1']))])))
        ig = ig_['O'+str(int(v_['3']))+','+str(int(v_['1']))]
        posel = len(self.grelt[ig])
        self.grelt  = (
              loaset(self.grelt,ig,posel,ie_['S'+str(int(v_['3']))+','+str(int(v_['1']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw  = (
              loaset(self.grelw,ig,posel,float(v_['T'+str(int(v_['3']))+','+str(int(v_['1']))])))
        posel = len(self.grelt[ig])
        self.grelt  = (
              loaset(self.grelt,ig,posel,ie_['V'+str(int(v_['3']))+','+str(int(v_['1']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw  = (
              loaset(self.grelw,ig,posel,float(v_['T'+str(int(v_['2']))+','+str(int(v_['1']))])))
        ig = ig_['O'+str(int(v_['2']))+','+str(int(v_['2']))]
        posel = len(self.grelt[ig])
        self.grelt  = (
              loaset(self.grelt,ig,posel,ie_['S'+str(int(v_['2']))+','+str(int(v_['2']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw  = (
              loaset(self.grelw,ig,posel,float(v_['T'+str(int(v_['2']))+','+str(int(v_['2']))])))
        posel = len(self.grelt[ig])
        self.grelt  = (
              loaset(self.grelt,ig,posel,ie_['U'+str(int(v_['2']))+','+str(int(v_['2']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw  = (
              loaset(self.grelw,ig,posel,float(v_['T'+str(int(v_['2']))+','+str(int(v_['1']))])))
        posel = len(self.grelt[ig])
        self.grelt  = (
              loaset(self.grelt,ig,posel,ie_['V'+str(int(v_['2']))+','+str(int(v_['2']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw  = (
              loaset(self.grelw,ig,posel,float(v_['T'+str(int(v_['2']))+','+str(int(v_['1']))])))
        posel = len(self.grelt[ig])
        self.grelt  = (
              loaset(self.grelt,ig,posel,ie_['W'+str(int(v_['2']))+','+str(int(v_['2']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw  = (
              loaset(self.grelw,ig,posel,float(v_['T'+str(int(v_['1']))+','+str(int(v_['1']))])))
        ig = ig_['O'+str(int(v_['3']))+','+str(int(v_['2']))]
        posel = len(self.grelt[ig])
        self.grelt  = (
              loaset(self.grelt,ig,posel,ie_['S'+str(int(v_['3']))+','+str(int(v_['2']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw  = (
              loaset(self.grelw,ig,posel,float(v_['T'+str(int(v_['3']))+','+str(int(v_['2']))])))
        posel = len(self.grelt[ig])
        self.grelt  = (
              loaset(self.grelt,ig,posel,ie_['U'+str(int(v_['3']))+','+str(int(v_['2']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw  = (
              loaset(self.grelw,ig,posel,float(v_['T'+str(int(v_['3']))+','+str(int(v_['1']))])))
        posel = len(self.grelt[ig])
        self.grelt  = (
              loaset(self.grelt,ig,posel,ie_['V'+str(int(v_['3']))+','+str(int(v_['2']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw  = (
              loaset(self.grelw,ig,posel,float(v_['T'+str(int(v_['2']))+','+str(int(v_['2']))])))
        posel = len(self.grelt[ig])
        self.grelt  = (
              loaset(self.grelt,ig,posel,ie_['W'+str(int(v_['3']))+','+str(int(v_['2']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw  = (
              loaset(self.grelw,ig,posel,float(v_['T'+str(int(v_['2']))+','+str(int(v_['1']))])))
        ig = ig_['O'+str(int(v_['3']))+','+str(int(v_['3']))]
        posel = len(self.grelt[ig])
        self.grelt  = (
              loaset(self.grelt,ig,posel,ie_['S'+str(int(v_['3']))+','+str(int(v_['3']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw  = (
              loaset(self.grelw,ig,posel,float(v_['T'+str(int(v_['3']))+','+str(int(v_['3']))])))
        posel = len(self.grelt[ig])
        self.grelt  = (
              loaset(self.grelt,ig,posel,ie_['U'+str(int(v_['3']))+','+str(int(v_['3']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw  = (
              loaset(self.grelw,ig,posel,float(v_['T'+str(int(v_['3']))+','+str(int(v_['2']))])))
        posel = len(self.grelt[ig])
        self.grelt  = (
              loaset(self.grelt,ig,posel,ie_['V'+str(int(v_['3']))+','+str(int(v_['3']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw  = (
              loaset(self.grelw,ig,posel,float(v_['T'+str(int(v_['3']))+','+str(int(v_['2']))])))
        posel = len(self.grelt[ig])
        self.grelt  = (
              loaset(self.grelt,ig,posel,ie_['W'+str(int(v_['3']))+','+str(int(v_['3']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw  = (
              loaset(self.grelw,ig,posel,float(v_['T'+str(int(v_['2']))+','+str(int(v_['2']))])))
        for I in range(int(v_['4']),int(v_['N'])+1):
            v_['I-1'] = -1+I
            v_['I-2'] = -2+I
            ig = ig_['O'+str(I)+','+str(int(v_['I-2']))]
            posel = len(self.grelt[ig])
            self.grelt  = (
                  loaset(self.grelt,ig,posel,ie_['S'+str(I)+','+str(int(v_['I-2']))]))
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw  = (
                  loaset(self.grelw,ig,posel,float(v_['T'+str(I)+','+str(int(v_['I-2']))])))
            posel = len(self.grelt[ig])
            self.grelt  = (
                  loaset(self.grelt,ig,posel,ie_['V'+str(I)+','+str(int(v_['I-2']))]))
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw  = (
                  loaset(self.grelw,ig,posel,float(v_['T'+str(int(v_['I-1']))+','+str(int(v_['I-2']))])))
            ig = ig_['O'+str(I)+','+str(int(v_['I-1']))]
            posel = len(self.grelt[ig])
            self.grelt  = (
                  loaset(self.grelt,ig,posel,ie_['S'+str(I)+','+str(int(v_['I-1']))]))
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw  = (
                  loaset(self.grelw,ig,posel,float(v_['T'+str(I)+','+str(int(v_['I-1']))])))
            posel = len(self.grelt[ig])
            self.grelt  = (
                  loaset(self.grelt,ig,posel,ie_['U'+str(I)+','+str(int(v_['I-1']))]))
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw  = (
                  loaset(self.grelw,ig,posel,float(v_['T'+str(I)+','+str(int(v_['I-2']))])))
            posel = len(self.grelt[ig])
            self.grelt  = (
                  loaset(self.grelt,ig,posel,ie_['V'+str(I)+','+str(int(v_['I-1']))]))
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw  = (
                  loaset(self.grelw,ig,posel,float(v_['T'+str(int(v_['I-1']))+','+str(int(v_['I-1']))])))
            posel = len(self.grelt[ig])
            self.grelt  = (
                  loaset(self.grelt,ig,posel,ie_['W'+str(I)+','+str(int(v_['I-1']))]))
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw  = (
                  loaset(self.grelw,ig,posel,float(v_['T'+str(int(v_['I-1']))+','+str(int(v_['I-2']))])))
            ig = ig_['O'+str(I)+','+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['S'+str(I)+','+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['T'+str(I)+','+str(I)]))
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['U'+str(I)+','+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw  = (
                  loaset(self.grelw,ig,posel,float(v_['T'+str(I)+','+str(int(v_['I-1']))])))
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['V'+str(I)+','+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw  = (
                  loaset(self.grelw,ig,posel,float(v_['T'+str(I)+','+str(int(v_['I-1']))])))
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['W'+str(I)+','+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw  = (
                  loaset(self.grelw,ig,posel,float(v_['T'+str(int(v_['I-1']))+','+str(int(v_['I-1']))])))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN(10)           6.00000000
# LO SOLTN(100)          68.0000000
# LO SOLTN(500)          340.000000
# LO SOLTN(1000)         ???
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CNOR2-AN-V-V"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def en2PR(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0,0]*EV_[1,0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1,0]
            g_[1] = EV_[0,0]
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

