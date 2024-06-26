from s2mpjlib import *
class  BROYDN3D(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : BROYDN3D
#    *********
# 
#    Broyden tridiagonal system on nonlinear equations
# 
#    Source:  problem 30 in
#    J.J. More', B.S. Garbow and K.E. Hillstrom,
#    "Testing Unconstrained Optimization Software",
#    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
# 
#    See also Toint#17 and Buckley#78.
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "NQR2-AN-V-V"
# 
#    N is the number of variables (variable).
# 
#           Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER     original value
# IE N                   50             $-PARAMETER
# IE N                   100            $-PARAMETER
# IE N                   500            $-PARAMETER
# IE N                   1000           $-PARAMETER
# IE N                   5000           $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'BROYDN3D'

    def __init__(self, *args): 
        import numpy as np
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(10);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
# IE N                   10000          $-PARAMETER
        if nargin<2:
            v_['KAPPA1'] = float(2.0);  #  SIF file default value
        else:
            v_['KAPPA1'] = float(args[1])
        if nargin<3:
            v_['KAPPA2'] = float(1.0);  #  SIF file default value
        else:
            v_['KAPPA2'] = float(args[2])
        v_['1'] = 1
        v_['2'] = 2
        v_['-K2'] = -1.0*v_['KAPPA2']
        v_['N-1'] = -1+v_['N']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.A       = lil_matrix((1000000,1000000))
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames      = np.array([])
        self.cnames = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('E1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E1')
        iv = ix_['X'+str(int(v_['2']))]
        self.A[ig,iv] = float(-2.0)+self.A[ig,iv]
        for I in range(int(v_['2']),int(v_['N-1'])+1):
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            [ig,ig_,_] = s2mpj_ii('E'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'E'+str(I))
            iv = ix_['X'+str(int(v_['I-1']))]
            self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
            iv = ix_['X'+str(int(v_['I+1']))]
            self.A[ig,iv] = float(-2.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('E'+str(int(v_['N'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E'+str(int(v_['N'])))
        iv = ix_['X'+str(int(v_['N-1']))]
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
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        for I in range(int(v_['1']),int(v_['N'])+1):
            self.gconst = arrset(self.gconst,ig_['E'+str(I)],float(v_['-K2']))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        self.xlower = np.zeros((self.n,1))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(-1.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eBROY', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftp = []
        elftp = loaset(elftp,it,0,'K1')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['1']),int(v_['N'])+1):
            ename = 'B'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eBROY')
            ielftype = arrset(ielftype, ie, iet_["eBROY"])
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,-1.0)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='K1')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['KAPPA1']))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            ig = ig_['E'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['B'+str(I)])
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
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        self.A.resize(ngrp,self.n)
        self.A     = self.A.tocsr()
        sA1,sA2    = self.A.shape
        self.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons =  np.where(self.congrps in np.setdiff1d(nlc,self.congrps))[0]
        self.pbclass = "NQR2-AN-V-V"
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eBROY(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        TEMP1 = -self.elpar[iel_][0]*EV_[0]
        TEMP2 = 3.0+TEMP1
        f_   = TEMP2*EV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = TEMP2+TEMP1
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = -self.elpar[iel_][0]-self.elpar[iel_][0]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

