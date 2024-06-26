from s2mpjlib import *
class  STANCMIN(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    A simple fractional programming problem.
# 
#    Source:
#    I.M. Stancu-Minasian
#    Private communication, 1992.
# 
#    SDIF input: Ph. Toint, October 1992.
# 
#    classification = "OLI2-AY-3-2"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'STANCMIN'

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
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.A       = lil_matrix((1000000,1000000))
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames      = np.array([])
        self.cnames = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('C1',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'C1')
        iv = ix_['X1']
        self.A[ig,iv] = float(3.0)+self.A[ig,iv]
        iv = ix_['X2']
        self.A[ig,iv] = float(4.0)+self.A[ig,iv]
        iv = ix_['X3']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C2',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'C2')
        iv = ix_['X1']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X2']
        self.A[ig,iv] = float(4.0)+self.A[ig,iv]
        iv = ix_['X3']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
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
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        self.gconst = arrset(self.gconst,ig_['C1'],float(2.0))
        self.gconst = arrset(self.gconst,ig_['C2'],float(1.0))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(50.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eRATIO', iet_)
        elftv = loaset(elftv,it,0,'XA')
        elftv = loaset(elftv,it,1,'XB')
        elftv = loaset(elftv,it,2,'XC')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        ename = 'E'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eRATIO')
        ielftype = arrset(ielftype, ie, iet_["eRATIO"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,50.0)
        posev = np.where(elftv[ielftype[ie]]=='XA')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,50.0)
        posev = np.where(elftv[ielftype[ie]]=='XB')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,50.0)
        posev = np.where(elftv[ielftype[ie]]=='XC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        ig = ig_['OBJ']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 5.0
#    Solution
# LO SOLTN               4.25
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.cupper[np.arange(self.nle)] = np.zeros((self.nle,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        self.A.resize(ngrp,self.n)
        self.A     = self.A.tocsr()
        sA1,sA2    = self.A.shape
        self.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons =  np.where(self.congrps in np.setdiff1d(nlc,self.congrps))[0]
        self.pbclass = "OLI2-AY-3-2"
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eRATIO(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((2,3))
        IV_ = np.zeros(2)
        U_[0,1] = U_[0,1]+6
        U_[0,0] = U_[0,0]+3
        U_[0,2] = U_[0,2]+2
        U_[1,0] = U_[1,0]+1
        U_[1,1] = U_[1,1]+4
        U_[1,2] = U_[1,2]+1
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        f_   = -(IV_[0]-11.0)/(IV_[1]+1.0)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -1.0/(IV_[1]+1.0)
            g_[1] = (IV_[0]-11.0)/(IV_[1]+1.0)**2
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = 1.0/(IV_[1]+1.0)**2
                H_[1,0] = H_[0,1]
                H_[1,1] = -2.0*(IV_[0]-11.0)/(IV_[1]+1.0)**3
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

