from s2mpjlib import *
class  LOGHAIRY(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : LOGHAIRY
#    *********
# 
#    A more difficult variant of the HAIRY problem in two variables.  
#    It is defined by a logarithmic transformation of the HAIRY surface,
#    which is defined by this function has a large number of relatively
#    sharp hills between  which a valley leads to the minimizer.
#    This problem contains a large number of saddle points.
# 
#    The problem is nonconvex.
# 
#    Source:
#    Ph. Toint, private communication,
# 
#    SIF input: Ph. Toint, April 1997.
# 
#    classification = "C-COUR2-AN-2-0"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 25 XI 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'LOGHAIRY'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['HLENGTH'] = 30.0
        v_['CSLOPE'] = 100.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [iv,ix_,_] = s2mpj_ii('X1',ix_)
        self.xnames=arrset(self.xnames,iv,'X1')
        [iv,ix_,_] = s2mpj_ii('X2',ix_)
        self.xnames=arrset(self.xnames,iv,'X2')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('FURCUP',ig_)
        gtype = arrset(gtype,ig,'<>')
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.x0[ix_['X1']] = float(-500.0)
        self.x0[ix_['X2']] = float(-700.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eFUR', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftp = []
        elftp = loaset(elftp,it,0,'DENS')
        [it,iet_,_] = s2mpj_ii( 'eDCUP', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftp = loaset(elftp,it,0,'SMOOTH')
        [it,iet_,_] = s2mpj_ii( 'en1CUP', iet_)
        elftv = loaset(elftv,it,0,'V')
        elftp = loaset(elftp,it,0,'SMOOTH')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        ename = 'HAIR'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eFUR')
        ielftype = arrset(ielftype,ie,iet_["eFUR"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='DENS')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(7.0))
        ename = 'DBOWL'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDCUP')
        ielftype = arrset(ielftype,ie,iet_["eDCUP"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='SMOOTH')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(0.01))
        ename = '1BOWL'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en1CUP')
        ielftype = arrset(ielftype,ie,iet_["en1CUP"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='V')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='SMOOTH')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(0.01))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gLOG',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        ig = ig_['FURCUP']
        self.grftype = arrset(self.grftype,ig,'gLOG')
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['HAIR'])
        self.grelw = loaset(self.grelw,ig,posel,float(v_['HLENGTH']))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['DBOWL'])
        self.grelw = loaset(self.grelw,ig,posel,float(v_['CSLOPE']))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['1BOWL'])
        self.grelw = loaset(self.grelw,ig,posel,float(v_['CSLOPE']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               0.1823216
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass   = "C-COUR2-AN-2-0"
        self.objderlvl = 2

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eFUR(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        DV1 = self.elpar[iel_][0]*EV_[0]
        DV2 = self.elpar[iel_][0]*EV_[1]
        TDV1 = DV1+DV1
        TDV2 = DV2+DV2
        TDL2 = 2.0*self.elpar[iel_][0]*self.elpar[iel_][0]
        S1SQ = np.sin(DV1)**2
        C2SQ = np.cos(DV2)**2
        STDV1 = np.sin(TDV1)
        STDV2 = np.sin(TDV2)
        f_   = S1SQ*C2SQ
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = self.elpar[iel_][0]*STDV1*C2SQ
            g_[1] = -self.elpar[iel_][0]*S1SQ*STDV2
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = TDL2*np.cos(TDV1)*C2SQ
                H_[0,1] = -self.elpar[iel_][0]*self.elpar[iel_][0]*STDV1*STDV2
                H_[1,0] = H_[0,1]
                H_[1,1] = -TDL2*S1SQ*np.cos(TDV2)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eDCUP(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((1,2))
        IV_ = np.zeros(1)
        U_[0,0] = U_[0,0]+1
        U_[0,1] = U_[0,1]-1
        IV_[0] = U_[0:1,:].dot(EV_)
        VSQ = IV_[0]*IV_[0]
        ARG = self.elpar[iel_][0]+VSQ
        SQARG = np.sqrt(ARG)
        DEN = 1.0/SQARG
        f_   = SQARG
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = IV_[0]*DEN
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = (1.0-VSQ/ARG)*DEN
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def en1CUP(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        VSQ = EV_[0]*EV_[0]
        ARG = self.elpar[iel_][0]+VSQ
        SQARG = np.sqrt(ARG)
        DEN = 1.0/SQARG
        f_   = SQARG
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[0]*DEN
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = (1.0-VSQ/ARG)*DEN
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gLOG(self,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        S = 1.0e2
        f_= np.log((S+GVAR_)/S)
        if nargout>1:
            g_ = 1.0/(S+GVAR_)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = -1.0/(S+GVAR_)**2
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

