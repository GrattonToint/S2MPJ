from s2mpjlib import *
class  BT9(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : BT9
#    *********
# 
#    Source: problem 9 in
#    P.T. Boggs and J.W. Tolle,
#    "A strategy for global convergence in a sequential 
#     quadratic programming algorithm",
#    SINUM 26(3), pp. 600-623, 1989.
# 
# 
#    The problem as stated in the paper seems to contain a typo.
#    In order to make the problem bounded below and the second constraint
#    feasible at the proposed solution, the sign of x2 in the second constraint
#    has been set to - instead of +.
# 
#    The problem is not convex.
# 
#    SIF input: Ph. Toint, June 1993.
# 
#    classification = "C-CLOR2-AN-4-2"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 17 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'BT9'

    def __init__(self, *args): 
        import numpy as np
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2mpj_ii('X1',ix_)
        self.xnames=arrset(self.xnames,iv,'X1')
        [iv,ix_,_] = s2mpj_ii('X2',ix_)
        self.xnames=arrset(self.xnames,iv,'X2')
        [iv,ix_,_] = s2mpj_ii('X3',ix_)
        self.xnames=arrset(self.xnames,iv,'X3')
        [iv,ix_,_] = s2mpj_ii('X4',ix_)
        self.xnames=arrset(self.xnames,iv,'X4')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.A       = lil_matrix((1000000,1000000))
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames      = np.array([])
        self.cnames = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X1']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('CON1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'CON1')
        iv = ix_['X2']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('CON2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'CON2')
        iv = ix_['X2']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
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
        self.cnames= cnames[self.congrps]
        self.nob = ngrp-self.m
        self.objgrps = np.where(gtype=='<>')[0]
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(2.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'X')
        [it,iet_,_] = s2mpj_ii( 'eCB', iet_)
        elftv = loaset(elftv,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        ename = 'X3SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQ')
        ielftype = arrset(ielftype,ie,iet_["eSQ"])
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(2.0))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'X1SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQ')
        ielftype = arrset(ielftype,ie,iet_["eSQ"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(2.0))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'X4SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQ')
        ielftype = arrset(ielftype,ie,iet_["eSQ"])
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(2.0))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'X1CB'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCB')
        ielftype = arrset(ielftype,ie,iet_["eCB"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(2.0))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        ig = ig_['CON1']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['X1CB'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['X3SQ'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        ig = ig_['CON2']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['X1SQ'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['X4SQ'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN                -1.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        self.A.resize(ngrp,self.n)
        self.A     = self.A.tocsr()
        sA1,sA2    = self.A.shape
        self.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass = "C-CLOR2-AN-4-2"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSQ(self, nargout,*args):

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
            g_[0] = 2.0*EV_[0]
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
    def eCB(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[0]*EV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 3.0*EV_[0]*EV_[0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 6.0*EV_[0]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

