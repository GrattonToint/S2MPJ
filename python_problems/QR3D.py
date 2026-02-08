from s2mpjlib import *
class  QR3D(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : QR3D
#    *********
# 
#    Find the QR factorization of a tridiagonal matrix A.
#    The problem is formulated as a system of quadratic equations
#    whose unknowns are the elements of the orthogonal matrix Q and of
#    the upper triangular matrix R.  In this version of the problem,
#    the banded structure of R is not imposed as a constraint. See problem
#    QR3DBD for the case where this structure is explicitly used.
# 
#    The problem is non-convex.
# 
#    Source:
#    Ph. Toint, private communication.
# 
#    SIF input: Ph. Toint, Nov 1993
# 
#    classification = "C-CNQR2-AN-V-V"
# 
#    Define the matrix order M  ( M >= 3 ).
#    There are M * ( 3M + 1) / 2 variables and equations.
# 
#           Alternative values for the SIF file parameters:
# IE M                   5              $-PARAMETER  n =  40
# IE M                   10             $-PARAMETER  n = 155  original value
# IE M                   20             $-PARAMETER  n = 610
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'QR3D'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['M'] = int(5);  #  SIF file default value
        else:
            v_['M'] = int(args[0])
        v_['1'] = 1
        v_['2'] = 2
        v_['M-1'] = -1+v_['M']
        v_['RM'] = float(v_['M'])
        v_['2/M'] = 2.0/v_['RM']
        v_['A'+str(int(v_['1']))+','+str(int(v_['1']))] = v_['2/M']
        v_['A'+str(int(v_['1']))+','+str(int(v_['2']))] = 0.0
        for I in range(int(v_['2']),int(v_['M-1'])+1):
            v_['I+1'] = 1+I
            v_['I-1'] = -1+I
            v_['1-I'] = -1*v_['I-1']
            v_['R1-I'] = float(v_['1-I'])
            v_['1-I/M'] = v_['R1-I']/v_['RM']
            v_['2I'] = 2*I
            v_['R2I'] = float(v_['2I'])
            v_['2I/M'] = v_['R2I']/v_['RM']
            v_['A'+str(I)+','+str(int(v_['I-1']))] = v_['1-I/M']
            v_['A'+str(I)+','+str(I)] = v_['2I/M']
            v_['A'+str(I)+','+str(int(v_['I+1']))] = v_['1-I/M']
        v_['RM-1'] = float(v_['M-1'])
        v_['1-M'] = -1.0*v_['RM-1']
        v_['1-M/M'] = v_['1-M']/v_['RM']
        v_['2M'] = 2.0*v_['RM']
        v_['A'+str(int(v_['M']))+','+str(int(v_['M-1']))] = v_['1-M/M']
        v_['A'+str(int(v_['M']))+','+str(int(v_['M']))] = v_['2M']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['M'])+1):
            for J in range(int(v_['1']),int(v_['M'])+1):
                [iv,ix_,_] = s2mpj_ii('Q'+str(I)+','+str(J),ix_)
                self.xnames=arrset(self.xnames,iv,'Q'+str(I)+','+str(J))
        for I in range(int(v_['1']),int(v_['M'])+1):
            for J in range(int(I),int(v_['M'])+1):
                [iv,ix_,_] = s2mpj_ii('R'+str(I)+','+str(J),ix_)
                self.xnames=arrset(self.xnames,iv,'R'+str(I)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['1']),int(v_['M'])+1):
            for J in range(int(I),int(v_['M'])+1):
                [ig,ig_,_] = s2mpj_ii('O'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'O'+str(I)+','+str(J))
        for I in range(int(v_['1']),int(v_['M'])+1):
            for J in range(int(v_['1']),int(v_['M'])+1):
                [ig,ig_,_] = s2mpj_ii('F'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'F'+str(I)+','+str(J))
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
        for I in range(int(v_['1']),int(v_['M'])+1):
            self.gconst = arrset(self.gconst,ig_['O'+str(I)+','+str(I)],float(1.0))
        self.gconst  = (
              arrset(self.gconst,ig_['F'+str(int(v_['1']))+','+str(int(v_['1']))],float(v_['A'+str(int(v_['1']))+','+str(int(v_['1']))])))
        self.gconst  = (
              arrset(self.gconst,ig_['F'+str(int(v_['1']))+','+str(int(v_['2']))],float(v_['A'+str(int(v_['1']))+','+str(int(v_['2']))])))
        for I in range(int(v_['2']),int(v_['M-1'])+1):
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            self.gconst  = (
                  arrset(self.gconst,ig_['F'+str(I)+','+str(int(v_['I-1']))],float(v_['A'+str(I)+','+str(int(v_['I-1']))])))
            self.gconst  = (
                  arrset(self.gconst,ig_['F'+str(I)+','+str(I)],float(v_['A'+str(I)+','+str(I)])))
            self.gconst  = (
                  arrset(self.gconst,ig_['F'+str(I)+','+str(int(v_['I+1']))],float(v_['A'+str(I)+','+str(int(v_['I+1']))])))
        self.gconst  = (
              arrset(self.gconst,ig_['F'+str(int(v_['M']))+','+str(int(v_['M-1']))],float(v_['A'+str(int(v_['M']))+','+str(int(v_['M-1']))])))
        self.gconst  = (
              arrset(self.gconst,ig_['F'+str(int(v_['M']))+','+str(int(v_['M']))],float(v_['A'+str(int(v_['M']))+','+str(int(v_['M']))])))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        for I in range(int(v_['1']),int(v_['M'])+1):
            self.xlower[ix_['R'+str(I)+','+str(I)]] = 0.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        for I in range(int(v_['1']),int(v_['M'])+1):
            self.x0[ix_['Q'+str(I)+','+str(I)]] = float(1.0)
        for I in range(int(v_['1']),int(v_['M-1'])+1):
            v_['I+1'] = 1+I
            self.x0[ix_['R'+str(I)+','+str(I)]] = float(v_['A'+str(I)+','+str(I)])
            self.x0[ix_['R'+str(I)+','+str(int(v_['I+1']))]]  = (
                  float(v_['A'+str(I)+','+str(int(v_['I+1']))]))
        self.x0[ix_['R'+str(int(v_['M']))+','+str(int(v_['M']))]]  = (
              float(v_['A'+str(int(v_['M']))+','+str(int(v_['M']))]))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['1']),int(v_['M'])+1):
            for J in range(int(I),int(v_['M'])+1):
                for K in range(int(v_['1']),int(v_['M'])+1):
                    ename = 'C'+str(I)+','+str(J)+','+str(K)
                    [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
                    if newelt:
                        self.elftype = arrset(self.elftype,ie,'en2PR')
                        ielftype = arrset(ielftype,ie,iet_['en2PR'])
                    vname = 'Q'+str(I)+','+str(K)
                    [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                    posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                    self.elvar = loaset(self.elvar,ie,posev[0],iv)
                    vname = 'Q'+str(J)+','+str(K)
                    [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                    posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                    self.elvar = loaset(self.elvar,ie,posev[0],iv)
        for I in range(int(v_['1']),int(v_['M'])+1):
            for J in range(int(v_['1']),int(v_['M'])+1):
                for K in range(int(v_['1']),int(J)+1):
                    ename = 'B'+str(I)+','+str(J)+','+str(K)
                    [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
                    if newelt:
                        self.elftype = arrset(self.elftype,ie,'en2PR')
                        ielftype = arrset(ielftype,ie,iet_['en2PR'])
                    vname = 'Q'+str(I)+','+str(K)
                    [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                    posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                    self.elvar = loaset(self.elvar,ie,posev[0],iv)
                    vname = 'R'+str(K)+','+str(J)
                    [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                    posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                    self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['M'])+1):
            for J in range(int(I),int(v_['M'])+1):
                for K in range(int(v_['1']),int(v_['M'])+1):
                    ig = ig_['O'+str(I)+','+str(J)]
                    posel = len(self.grelt[ig])
                    self.grelt  = (
                          loaset(self.grelt,ig,posel,ie_['C'+str(I)+','+str(J)+','+str(K)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    self.grelw = loaset(self.grelw,ig,posel,1.)
        for I in range(int(v_['1']),int(v_['M'])+1):
            for J in range(int(v_['1']),int(v_['M'])+1):
                for K in range(int(v_['1']),int(J)+1):
                    ig = ig_['F'+str(I)+','+str(J)]
                    posel = len(self.grelt[ig])
                    self.grelt  = (
                          loaset(self.grelt,ig,posel,ie_['B'+str(I)+','+str(J)+','+str(K)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CNQR2-AN-V-V"
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

