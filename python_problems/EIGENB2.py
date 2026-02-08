from s2mpjlib import *
class  EIGENB2(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : EIGENB2
#    --------
# 
#    Solving symmetric eigenvalue problems as systems of
#    nonlinear equations.
# 
#    The problem is, given a symmetric matrix A, to find an orthogonal
#    matrix Q and diagonal matrix D such that A Q(T) = Q(T) D.
# 
#    Example B: a tridiagonal matrix with diagonals 2 and off diagonals -1
# 
#    Source:  An idea by Nick Gould
# 
#                Nonlinear equations version 2.
# 
#    SIF input: Nick Gould, Nov 1992.
# 
#    classification = "C-CQQR2-AN-V-V"
# 
#    The dimension of the matrix.
# 
#           Alternative values for the SIF file parameters:
# IE N                   2              $-PARAMETER
# IE N                   10             $-PARAMETER     original value
# IE N                   50             $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'EIGENB2'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(2);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
# IE N                   100            $-PARAMETER
        v_['1'] = 1
        for J in range(int(v_['1']),int(v_['N'])+1):
            for I in range(int(v_['1']),int(v_['N'])+1):
                v_['A'+str(I)+','+str(J)] = 0.0
            v_['J-1'] = -1+J
            v_['A'+str(int(v_['J-1']))+','+str(J)] = -1.0
            v_['A'+str(J)+','+str(J)] = 2.0
            v_['J+1'] = 1+J
            v_['A'+str(int(v_['J+1']))+','+str(J)] = -1.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for J in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('D'+str(J),ix_)
            self.xnames=arrset(self.xnames,iv,'D'+str(J))
            for I in range(int(v_['1']),int(v_['N'])+1):
                [iv,ix_,_] = s2mpj_ii('Q'+str(I)+','+str(J),ix_)
                self.xnames=arrset(self.xnames,iv,'Q'+str(I)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for J in range(int(v_['1']),int(v_['N'])+1):
            for I in range(int(v_['1']),int(J)+1):
                [ig,ig_,_] = s2mpj_ii('O'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'O'+str(I)+','+str(J))
            for I in range(int(v_['1']),int(v_['N'])+1):
                for K in range(int(v_['1']),int(v_['N'])+1):
                    v_['-AIK'] = -1.0*v_['A'+str(I)+','+str(K)]
                    [ig,ig_,_] = s2mpj_ii('E'+str(I)+','+str(J),ig_)
                    gtype = arrset(gtype,ig,'<>')
                    irA  = np.append(irA,[ig])
                    icA  = np.append(icA,[ix_['Q'+str(J)+','+str(K)]])
                    valA = np.append(valA,float(v_['-AIK']))
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
        for J in range(int(v_['1']),int(v_['N'])+1):
            self.gconst = arrset(self.gconst,ig_['O'+str(J)+','+str(J)],float(1.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(0.0))
        for J in range(int(v_['1']),int(v_['N'])+1):
            v_['RJ'] = float(J)
            self.x0[ix_['D'+str(J)]] = float(1.0)
            self.x0[ix_['Q'+str(J)+','+str(J)]] = float(1.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'en2PROD', iet_)
        elftv = loaset(elftv,it,0,'Q1')
        elftv = loaset(elftv,it,1,'Q2')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for J in range(int(v_['1']),int(v_['N'])+1):
            for I in range(int(v_['1']),int(v_['N'])+1):
                ename = 'E'+str(I)+','+str(J)
                [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
                if newelt:
                    self.elftype = arrset(self.elftype,ie,'en2PROD')
                    ielftype = arrset(ielftype,ie,iet_['en2PROD'])
                vname = 'Q'+str(J)+','+str(I)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
                posev = np.where(elftv[ielftype[ie]]=='Q1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'D'+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
                posev = np.where(elftv[ielftype[ie]]=='Q2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
            for I in range(int(v_['1']),int(J)+1):
                for K in range(int(v_['1']),int(v_['N'])+1):
                    ename = 'O'+str(I)+','+str(J)+','+str(K)
                    [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
                    if newelt:
                        self.elftype = arrset(self.elftype,ie,'en2PROD')
                        ielftype = arrset(ielftype,ie,iet_['en2PROD'])
                    vname = 'Q'+str(K)+','+str(I)
                    [iv,ix_]  = (
                          s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0)))
                    posev = np.where(elftv[ielftype[ie]]=='Q1')[0]
                    self.elvar = loaset(self.elvar,ie,posev[0],iv)
                    vname = 'Q'+str(K)+','+str(J)
                    [iv,ix_]  = (
                          s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0)))
                    posev = np.where(elftv[ielftype[ie]]=='Q2')[0]
                    self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gL2',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for J in range(int(v_['1']),int(v_['N'])+1):
            for I in range(int(v_['1']),int(v_['N'])+1):
                ig = ig_['E'+str(I)+','+str(J)]
                self.grftype = arrset(self.grftype,ig,'gL2')
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,1.)
            for I in range(int(v_['1']),int(J)+1):
                for K in range(int(v_['1']),int(v_['N'])+1):
                    ig = ig_['O'+str(I)+','+str(J)]
                    posel = len(self.grelt[ig])
                    self.grelt  = (
                          loaset(self.grelt,ig,posel,ie_['O'+str(I)+','+str(J)+','+str(K)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CQQR2-AN-V-V"
        self.objderlvl = 2
        self.conderlvl = [2]


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def en2PROD(self, nargout,*args):

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
                H_[0,1] = 1.0e+0
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gL2(self,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= GVAR_*GVAR_
        if nargout>1:
            g_ = GVAR_+GVAR_
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 2.0e+0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

