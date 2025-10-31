from s2mpjlib import *
class  HADAMALS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HADAMALS
#    --------
# 
#    An attempt to find Hadamard matrices of order N.
# 
#    The problem is to find an N by N orthonormal matrix Q,
#    with column norms N, whose entries are plus or minus one.
# 
#    Source:  A suggestion by Alan Edelman (MIT).
# 
#    SIF input: Nick Gould, Nov 1993.
# 
#    classification = "C-COBR2-RN-V-V"
# 
#    The dimension of the matrix (=> N**2 variables).
# 
#           Alternative values for the SIF file parameters:
# IE N                   2              $-PARAMETER    original value
# IE N                   4              $-PARAMETER
# IE N                   6              $-PARAMETER
# IE N                   8              $-PARAMETER
# IE N                   10             $-PARAMETER
# IE N                   12             $-PARAMETER
# IE N                   14             $-PARAMETER
# IE N                   16             $-PARAMETER
# IE N                   18             $-PARAMETER
# IE N                   20             $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HADAMALS'

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
# IE N                   32             $-PARAMETER
# IE N                   64             $-PARAMETER
# IE N                   128            $-PARAMETER
# IE N                   256            $-PARAMETER
# IE N                   428            $-PARAMETER
        v_['1'] = 1
        v_['2'] = 2
        v_['RN'] = float(v_['N'])
        v_['N/2'] = int(np.fix(v_['N']/v_['2']))
        v_['N/2+1'] = 1+v_['N/2']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for J in range(int(v_['1']),int(v_['N'])+1):
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
                gtype = arrset(gtype,ig,'<>')
        for J in range(int(v_['1']),int(v_['N'])+1):
            for I in range(int(v_['2']),int(v_['N'])+1):
                [ig,ig_,_] = s2mpj_ii('S'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        for J in range(int(v_['1']),int(v_['N'])+1):
            self.gconst = arrset(self.gconst,ig_['O'+str(J)+','+str(J)],float(v_['RN']))
        for J in range(int(v_['1']),int(v_['N'])+1):
            for I in range(int(v_['2']),int(v_['N'])+1):
                self.gconst = arrset(self.gconst,ig_['S'+str(I)+','+str(J)],float(1.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-1.0)
        self.xupper = np.full((self.n,1),1.0)
        for I in range(int(v_['1']),int(v_['N/2'])+1):
            self.xlower[ix_['Q'+str(I)+','+str(int(v_['1']))]] = 1.0
            self.xupper[ix_['Q'+str(I)+','+str(int(v_['1']))]] = 1.0
        for I in range(int(v_['N/2+1']),int(v_['N'])+1):
            self.xlower[ix_['Q'+str(I)+','+str(int(v_['1']))]] = -1.0
            self.xupper[ix_['Q'+str(I)+','+str(int(v_['1']))]] = -1.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        for J in range(int(v_['1']),int(v_['N'])+1):
            for I in range(int(v_['1']),int(v_['N/2'])+1):
                self.x0[ix_['Q'+str(I)+','+str(J)]] = float(0.9)
            for I in range(int(v_['N/2+1']),int(v_['N'])+1):
                self.x0[ix_['Q'+str(I)+','+str(J)]] = float(-0.9)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQR', iet_)
        elftv = loaset(elftv,it,0,'Q1')
        [it,iet_,_] = s2mpj_ii( 'en2PROD', iet_)
        elftv = loaset(elftv,it,0,'Q1')
        elftv = loaset(elftv,it,1,'Q2')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for J in range(int(v_['1']),int(v_['N'])+1):
            for I in range(int(v_['1']),int(J)+1):
                for K in range(int(v_['1']),int(v_['N'])+1):
                    ename = 'O'+str(I)+','+str(J)+','+str(K)
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    self.elftype = arrset(self.elftype,ie,'en2PROD')
                    ielftype = arrset(ielftype,ie,iet_["en2PROD"])
                    vname = 'Q'+str(K)+','+str(I)
                    [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-1.0),float(1.0),None)
                    posev = np.where(elftv[ielftype[ie]]=='Q1')[0]
                    self.elvar = loaset(self.elvar,ie,posev[0],iv)
                    vname = 'Q'+str(K)+','+str(J)
                    [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-1.0),float(1.0),None)
                    posev = np.where(elftv[ielftype[ie]]=='Q2')[0]
                    self.elvar = loaset(self.elvar,ie,posev[0],iv)
        for J in range(int(v_['1']),int(v_['N'])+1):
            for I in range(int(v_['2']),int(v_['N'])+1):
                ename = 'S'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eSQR')
                ielftype = arrset(ielftype,ie,iet_["eSQR"])
                vname = 'Q'+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-1.0),float(1.0),None)
                posev = np.where(elftv[ielftype[ie]]=='Q1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gL2',igt_)
        [it,igt_,_] = s2mpj_ii('gLARGEL2',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for J in range(int(v_['1']),int(v_['N'])+1):
            for I in range(int(v_['1']),int(J)+1):
                ig = ig_['O'+str(I)+','+str(J)]
                self.grftype = arrset(self.grftype,ig,'gL2')
                for K in range(int(v_['1']),int(v_['N'])+1):
                    ig = ig_['O'+str(I)+','+str(J)]
                    posel = len(self.grelt[ig])
                    self.grelt  = (
                          loaset(self.grelt,ig,posel,ie_['O'+str(I)+','+str(J)+','+str(K)]))
                    self.grelw = loaset(self.grelw,ig,posel,1.)
        for J in range(int(v_['1']),int(v_['N'])+1):
            for I in range(int(v_['2']),int(v_['N'])+1):
                ig = ig_['S'+str(I)+','+str(J)]
                self.grftype = arrset(self.grftype,ig,'gLARGEL2')
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['S'+str(I)+','+str(J)])
                self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass   = "C-COBR2-RN-V-V"
        self.objderlvl = 2


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSQR(self, nargout,*args):

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
                H_[0,0] = 2.0e+0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def en2PROD(self, nargout,*args):

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

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def g_globs(self):

        self.gfpar = np.array([]);
        self.gfpar = arrset(self.gfpar,0,1.0e+0)     # this is FACTOR
        return pbm

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

    @staticmethod
    def gLARGEL2(self,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= self.gfpar[0]*GVAR_*GVAR_
        if nargout>1:
            g_ = 2.0e+0*self.gfpar[0]*GVAR_
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 2.0e+0*self.gfpar[0]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

