from s2mpjlib import *
class  SPIRAL(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : SPIRAL
#    *********
# 
#    A nonlinear minmax problem.
# 
#    Source:
#    E. Polak, J.E. Higgins and D. Mayne,
#    "A barrier function for minmax problems",
#    Mathematical Programming, vol.54(2), pp. 155-176, 1992.
# 
#    SIF input: Ph. Toint, April 1992.
# 
#    classification = "C-CLOR2-AN-3-2"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'SPIRAL'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
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
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [iv,ix_,_] = s2mpj_ii('X1',ix_)
        self.xnames=arrset(self.xnames,iv,'X1')
        [iv,ix_,_] = s2mpj_ii('X2',ix_)
        self.xnames=arrset(self.xnames,iv,'X2')
        [iv,ix_,_] = s2mpj_ii('U',ix_)
        self.xnames=arrset(self.xnames,iv,'U')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['U']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('C1',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C1')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['U']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('C2',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C2')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['U']])
        valA = np.append(valA,float(1.0))
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
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        self.x0[ix_['X1']] = float(1.41831)
        self.x0[ix_['X2']] = float(-4.79462)
        self.x0[ix_['U']] = float(1.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'X')
        [it,iet_,_] = s2mpj_ii( 'eBADCOS', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        [it,iet_,_] = s2mpj_ii( 'eBADSIN', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        ename = 'X1SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQ')
        ielftype = arrset(ielftype,ie,iet_["eSQ"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'X2SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQ')
        ielftype = arrset(ielftype,ie,iet_["eSQ"])
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'BC'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eBADCOS')
        ielftype = arrset(ielftype,ie,iet_["eBADCOS"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'BS'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eBADSIN')
        ielftype = arrset(ielftype,ie,iet_["eBADSIN"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        ig = ig_['C1']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['BC'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['X1SQ'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-0.005))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['X2SQ'])
        self.grelw = loaset(self.grelw,ig,posel,float(-0.005))
        ig = ig_['C2']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['BS'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['X1SQ'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-0.005))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['X2SQ'])
        self.grelw = loaset(self.grelw,ig,posel,float(-0.005))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN              0.0
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle+self.neq,self.m)] = np.zeros((self.nge,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CLOR2-AN-3-2"
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
        f_   = EV_[0,0]*EV_[0,0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[0,0]+EV_[0,0]
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
    def eBADCOS(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        XX = EV_[0,0]*EV_[0,0]
        YY = EV_[1,0]*EV_[1,0]
        R = np.sqrt(XX+YY)
        DRDX = EV_[0,0]/R
        DRDY = EV_[1,0]/R
        R3 = R**3
        D2RDXX = 1.0/R-XX/R3
        D2RDYY = 1.0/R-YY/R3
        D2RDXY = -EV_[0,0]*EV_[1,0]/R3
        C = np.cos(R)
        S = np.sin(R)
        DCDX = -S*DRDX
        DCDY = -S*DRDY
        D2CDXX = -C*DRDX*DRDX-S*D2RDXX
        D2CDYY = -C*DRDY*DRDY-S*D2RDYY
        D2CDXY = -C*DRDX*DRDY-S*D2RDXY
        DSDX = C*DRDX
        DSDY = C*DRDY
        D2SDXX = -S*DRDX*DRDX+C*D2RDXX
        D2SDYY = -S*DRDY*DRDY+C*D2RDYY
        D2SDXY = -S*DRDX*DRDY+C*D2RDXY
        Z = EV_[0,0]-R*C
        DZDX = 1.0-DRDX*C-R*DCDX
        DZDY = -DRDY*C-R*DCDY
        D2ZDXX = -D2RDXX*C-2.0*DRDX*DCDX-R*D2CDXX
        D2ZDYY = -D2RDYY*C-2.0*DRDY*DCDY-R*D2CDYY
        D2ZDXY = -D2RDXY*C-DRDX*DCDY-DRDY*DCDX-R*D2CDXY
        f_   = Z*Z
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*DZDX*Z
            g_[1] = 2.0*DZDY*Z
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 2.0*(D2ZDXX*Z+DZDX*DZDX)
                H_[0,1] = 2.0*(D2ZDXY*Z+DZDX*DZDY)
                H_[1,0] = H_[0,1]
                H_[1,1] = 2.0*(D2ZDYY*Z+DZDY*DZDY)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eBADSIN(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        XX = EV_[0,0]*EV_[0,0]
        YY = EV_[1,0]*EV_[1,0]
        R = np.sqrt(XX+YY)
        DRDX = EV_[0,0]/R
        DRDY = EV_[1,0]/R
        R3 = R**3
        D2RDXX = 1.0/R-XX/R3
        D2RDYY = 1.0/R-YY/R3
        D2RDXY = -EV_[0,0]*EV_[1,0]/R3
        C = np.cos(R)
        S = np.sin(R)
        DCDX = -S*DRDX
        DCDY = -S*DRDY
        D2CDXX = -C*DRDX*DRDX-S*D2RDXX
        D2CDYY = -C*DRDY*DRDY-S*D2RDYY
        D2CDXY = -C*DRDX*DRDY-S*D2RDXY
        DSDX = C*DRDX
        DSDY = C*DRDY
        D2SDXX = -S*DRDX*DRDX+C*D2RDXX
        D2SDYY = -S*DRDY*DRDY+C*D2RDYY
        D2SDXY = -S*DRDX*DRDY+C*D2RDXY
        Z = EV_[1,0]-R*S
        DZDX = -DRDX*S-R*DSDX
        DZDY = 1.0-DRDY*S-R*DSDY
        D2ZDXX = -D2RDXX*S-2.0*DRDX*DSDX-R*D2SDXX
        D2ZDYY = -D2RDYY*S-2.0*DRDY*DSDY-R*D2SDYY
        D2ZDXY = -D2RDXY*S-DRDX*DSDY-DRDY*DSDX-R*D2SDXY
        f_   = Z*Z
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*DZDX*Z
            g_[1] = 2.0*DZDY*Z
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 2.0*(D2ZDXX*Z+DZDX*DZDX)
                H_[0,1] = 2.0*(D2ZDXY*Z+DZDX*DZDY)
                H_[1,0] = H_[0,1]
                H_[1,1] = 2.0*(D2ZDYY*Z+DZDY*DZDY)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

