from s2mpjlib import *
class  SYNTHES2(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : SYNTHES2
#    *********
# 
#    Source: Test problem 2 (Synthesis of processing system) in 
#    M. Duran & I.E. Grossmann,
#    "An outer approximation algorithm for a class of mixed integer nonlinear
#     programs", Mathematical Programming 36, pp. 307-339, 1986.
# 
#    SIF input: S. Leyffer, October 1997
# 
#    classification = "C-COOR2-AN-11-14"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 17 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'SYNTHES2'

    def __init__(self, *args): 
        import numpy as np
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['1'] = 1
        v_['5'] = 5
        v_['6'] = 6
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['6'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I))
        for I in range(int(v_['1']),int(v_['5'])+1):
            [iv,ix_,_] = s2mpj_ii('Y'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'Y'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.A       = lil_matrix((1000000,1000000))
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames      = np.array([])
        self.cnames = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['Y1']
        self.A[ig,iv] = float(5.0)+self.A[ig,iv]
        iv = ix_['Y2']
        self.A[ig,iv] = float(8.0)+self.A[ig,iv]
        iv = ix_['Y3']
        self.A[ig,iv] = float(6.0)+self.A[ig,iv]
        iv = ix_['Y4']
        self.A[ig,iv] = float(10.0)+self.A[ig,iv]
        iv = ix_['Y5']
        self.A[ig,iv] = float(6.0)+self.A[ig,iv]
        iv = ix_['X1']
        self.A[ig,iv] = float(-10.0)+self.A[ig,iv]
        iv = ix_['X2']
        self.A[ig,iv] = float(-15.0)+self.A[ig,iv]
        iv = ix_['X3']
        self.A[ig,iv] = float(-15.0)+self.A[ig,iv]
        iv = ix_['X4']
        self.A[ig,iv] = float(15.0)+self.A[ig,iv]
        iv = ix_['X5']
        self.A[ig,iv] = float(5.0)+self.A[ig,iv]
        iv = ix_['X6']
        self.A[ig,iv] = float(-20.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('N1',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'N1')
        [ig,ig_,_] = s2mpj_ii('N2',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'N2')
        iv = ix_['Y1']
        self.A[ig,iv] = float(-10.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('N3',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'N3')
        iv = ix_['Y2']
        self.A[ig,iv] = float(-10.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('L1',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'L1')
        iv = ix_['X3']
        self.A[ig,iv] = float(1.25)+self.A[ig,iv]
        iv = ix_['Y3']
        self.A[ig,iv] = float(-10.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('L2',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'L2')
        iv = ix_['X4']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X5']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['Y4']
        self.A[ig,iv] = float(-10.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('L3',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'L3')
        iv = ix_['X3']
        self.A[ig,iv] = float(-2.0)+self.A[ig,iv]
        iv = ix_['X6']
        self.A[ig,iv] = float(2.0)+self.A[ig,iv]
        iv = ix_['Y5']
        self.A[ig,iv] = float(-10.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('L4',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'L4')
        iv = ix_['X1']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X2']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X3']
        self.A[ig,iv] = float(-2.0)+self.A[ig,iv]
        iv = ix_['X4']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X6']
        self.A[ig,iv] = float(2.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('L5',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'L5')
        iv = ix_['X1']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X2']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X3']
        self.A[ig,iv] = float(-0.75)+self.A[ig,iv]
        iv = ix_['X4']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X6']
        self.A[ig,iv] = float(2.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('L6',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'L6')
        iv = ix_['X3']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X6']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('L7',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'L7')
        iv = ix_['X3']
        self.A[ig,iv] = float(2.0)+self.A[ig,iv]
        iv = ix_['X4']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X6']
        self.A[ig,iv] = float(-2.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('L8',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'L8')
        iv = ix_['X4']
        self.A[ig,iv] = float(-0.5)+self.A[ig,iv]
        iv = ix_['X5']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('L9',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'L9')
        iv = ix_['X4']
        self.A[ig,iv] = float(0.2)+self.A[ig,iv]
        iv = ix_['X5']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('L10',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'L10')
        iv = ix_['Y1']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['Y2']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('L11',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'L11')
        iv = ix_['Y4']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['Y5']
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
        self.gconst = arrset(self.gconst,ig_['OBJ'],float(-140.0))
        self.gconst = arrset(self.gconst,ig_['N2'],float(1.0))
        self.gconst = arrset(self.gconst,ig_['N3'],float(1.0))
        self.gconst = arrset(self.gconst,ig_['L10'],float(1.0))
        self.gconst = arrset(self.gconst,ig_['L11'],float(1.0))
        #%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange = np.full((ngrp,1),None)
        grange[legrps] = np.full((self.nle,1),float('inf'))
        grange = arrset(grange,ig_['L11'],float(-1.0))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        self.xupper[ix_['X1']] = 2.0
        self.xupper[ix_['X2']] = 2.0
        self.xupper[ix_['X3']] = 2.0
        self.xupper[ix_['X6']] = 3.0
        for I in range(int(v_['1']),int(v_['5'])+1):
            self.xupper[ix_['Y'+str(I)]] = 1.0
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eLOGSUM', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        [it,iet_,_] = s2mpj_ii( 'eEXPA', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftp = []
        elftp = loaset(elftp,it,0,'A')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        ename = 'LOGX4X5'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eLOGSUM')
        ielftype = arrset(ielftype,ie,iet_["eLOGSUM"])
        self.x0 = np.zeros((self.n,1))
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'EXPX1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eEXPA')
        ielftype = arrset(ielftype,ie,iet_["eEXPA"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='A')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
        ename = 'EXPX2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eEXPA')
        ielftype = arrset(ielftype,ie,iet_["eEXPA"])
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='A')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.2))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        ig = ig_['OBJ']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['EXPX1'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['EXPX2'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['LOGX4X5'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-60.0))
        ig = ig_['N1']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['LOGX4X5'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        ig = ig_['N2']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['EXPX1'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        ig = ig_['N3']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['EXPX2'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle)] = grange[legrps]
        self.cupper[np.arange(self.nle)] = np.zeros((self.nle,1))
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
        self.pbclass = "C-COOR2-AN-11-14"
        self.x0        = np.zeros((self.n,1))
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eLOGSUM(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = np.log(EV_[0]+EV_[1]+1.0)
        if not isinstance( f_, float ):
            f_   = f_.item();
        DX = 1.0/(EV_[0]+EV_[1]+1.0)
        DXDX = -1.0/(EV_[0]+EV_[1]+1.0)**2
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = DX
            g_[1] = DX
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = DXDX
                H_[0,1] = DXDX
                H_[1,0] = H_[0,1]
                H_[1,1] = DXDX
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eEXPA(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        EXPXA = np.exp(EV_[0]/self.elpar[iel_][0])
        f_   = EXPXA
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EXPXA/self.elpar[iel_][0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = EXPXA/self.elpar[iel_][0]/self.elpar[iel_][0]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

