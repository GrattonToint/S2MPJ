from s2mpjlib import *
class  LUKSAN14LS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : LUKSAN14LS
#    *********
# 
#    Problem 14 (chained and modified HS53) in the paper
# 
#      L. Luksan
#      Hybrid methods in large sparse nonlinear least squares
#      J. Optimization Theory & Applications 89(3) 575-595 (1996)
# 
#    SIF input: Nick Gould, June 2017.
# 
#    least-squares version
# 
#    classification = "SUR2-AN-V-0"
# 
#   seed for dimensions
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'LUKSAN14LS'

    def __init__(self, *args): 
        import numpy as np
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['S'] = 32
        v_['N'] = 3*v_['S']
        v_['N'] = 2+v_['N']
        v_['M'] = 7*v_['S']
        v_['1'] = 1
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
        v_['I'] = 1
        v_['K'] = 1
        for J in range(int(v_['1']),int(v_['S'])+1):
            v_['K+1'] = 1+v_['K']
            v_['K+2'] = 2+v_['K']
            v_['K+3'] = 3+v_['K']
            v_['K+4'] = 4+v_['K']
            v_['K+5'] = 5+v_['K']
            v_['K+6'] = 6+v_['K']
            v_['I+1'] = 1+v_['I']
            v_['I+2'] = 2+v_['I']
            v_['I+3'] = 3+v_['I']
            v_['I+4'] = 4+v_['I']
            [ig,ig_,_] = s2mpj_ii('E'+str(int(v_['K'])),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(int(v_['I+1']))]
            self.A[ig,iv] = float(-10.0e0)+self.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('E'+str(int(v_['K+1'])),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(int(v_['I+1']))]
            self.A[ig,iv] = float(1.0e0)+self.A[ig,iv]
            iv = ix_['X'+str(int(v_['I+2']))]
            self.A[ig,iv] = float(1.0e0)+self.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('E'+str(int(v_['K+2'])),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(int(v_['I+3']))]
            self.A[ig,iv] = float(1.0e0)+self.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('E'+str(int(v_['K+3'])),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(int(v_['I+4']))]
            self.A[ig,iv] = float(1.0e0)+self.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('E'+str(int(v_['K+4'])),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(int(v_['I']))]
            self.A[ig,iv] = float(1.0e0)+self.A[ig,iv]
            iv = ix_['X'+str(int(v_['I+1']))]
            self.A[ig,iv] = float(3.0e0)+self.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('E'+str(int(v_['K+5'])),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(int(v_['I+2']))]
            self.A[ig,iv] = float(1.0e0)+self.A[ig,iv]
            iv = ix_['X'+str(int(v_['I+3']))]
            self.A[ig,iv] = float(1.0e0)+self.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('E'+str(int(v_['K+5'])),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(int(v_['I+4']))]
            self.A[ig,iv] = float(-2.0e0)+self.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('E'+str(int(v_['K+6'])),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(int(v_['I+4']))]
            self.A[ig,iv] = float(-10.0e0)+self.A[ig,iv]
            v_['I'] = 3+v_['I']
            v_['K'] = 7+v_['K']
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        v_['K'] = 1
        for J in range(int(v_['1']),int(v_['S'])+1):
            v_['K+1'] = 1+v_['K']
            v_['K+2'] = 2+v_['K']
            v_['K+3'] = 3+v_['K']
            self.gconst = arrset(self.gconst,ig_['E'+str(int(v_['K+1']))],float(2.0))
            self.gconst = arrset(self.gconst,ig_['E'+str(int(v_['K+2']))],float(1.0))
            self.gconst = arrset(self.gconst,ig_['E'+str(int(v_['K+3']))],float(1.0))
            v_['K'] = 7+v_['K']
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        self.xlower = np.zeros((self.n,1))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        for I in range(int(v_['1']),int(v_['N'])+1):
            self.x0[ix_['X'+str(I)]] = float(-1.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQR', iet_)
        elftv = loaset(elftv,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        v_['I'] = 1
        v_['K'] = 1
        for J in range(int(v_['1']),int(v_['S'])+1):
            v_['K+6'] = 6+v_['K']
            v_['I+1'] = 1+v_['I']
            ename = 'E'+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype, ie, iet_["eSQR"])
            ename = 'E'+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(int(v_['I']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'E'+str(int(v_['K+6']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype, ie, iet_["eSQR"])
            ename = 'E'+str(int(v_['K+6']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(int(v_['I+1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            v_['I'] = 3+v_['I']
            v_['K'] = 7+v_['K']
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
        for ig in range(0,ngrp):
            self.grftype = arrset(self.grftype,ig,'gL2')
        v_['K'] = 1
        for J in range(int(v_['1']),int(v_['S'])+1):
            v_['K+6'] = 6+v_['K']
            ig = ig_['E'+str(int(v_['K']))]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['K']))])
            self.grelw = loaset(self.grelw,ig,posel,float(10.0))
            ig = ig_['E'+str(int(v_['K+6']))]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['K+6']))])
            self.grelw = loaset(self.grelw,ig,posel,float(10.0))
            v_['K'] = 7+v_['K']
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 0.0
#    Solution
# LO SOLTN                0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        self.A.resize(ngrp,self.n)
        self.A     = self.A.tocsr()
        sA1,sA2    = self.A.shape
        self.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass = "SUR2-AN-V-0"
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

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
                H_[0,0] = 2.0e0
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
                H_ = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

