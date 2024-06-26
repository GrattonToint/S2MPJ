from s2mpjlib import *
class  S268(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : S268
#    *********
# 
#    A quadratic programming problem.
# 
#    Source:
#    K. Schittkowski
#    "More Test Examples for Nonlinear Programming Codes"
#    Springer Verlag, Berlin, Lecture notes in economics and 
#    mathematical systems, volume 282, 1987
# 
#    SIF input: Michel Bierlaire and Annick Sartenaer,
#    October 1992.
# 
#    classification = "QLR2-AN-5-5"
# 
#   the number of functions
# 
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'S268'

    def __init__(self, *args): 
        import numpy as np
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['D1,1'] = 10197.0
        v_['D1,2'] = -12454.0
        v_['D1,3'] = -1013.0
        v_['D1,4'] = 1948.0
        v_['D1,5'] = 329.0
        v_['D2,1'] = -12454.0
        v_['D2,2'] = 20909.0
        v_['D2,3'] = -1733.0
        v_['D2,4'] = -4914.0
        v_['D2,5'] = -186.0
        v_['D3,1'] = -1013.0
        v_['D3,2'] = -1733.0
        v_['D3,3'] = 1755.0
        v_['D3,4'] = 1089.0
        v_['D3,5'] = -174.0
        v_['D4,1'] = 1948.0
        v_['D4,2'] = -4914.0
        v_['D4,3'] = 1089.0
        v_['D4,4'] = 1515.0
        v_['D4,5'] = -22.0
        v_['D5,1'] = 329.0
        v_['D5,2'] = -186.0
        v_['D5,3'] = -174.0
        v_['D5,4'] = -22.0
        v_['D5,5'] = 27.0
        v_['B1'] = -9170
        v_['B2'] = 17099
        v_['B3'] = -2271
        v_['B4'] = -4336
        v_['B5'] = -43
        v_['1'] = 1
        v_['5'] = 5
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['5'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.A       = lil_matrix((1000000,1000000))
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames      = np.array([])
        self.cnames = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('NONL',ig_)
        gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['1']),int(v_['5'])+1):
            [ig,ig_,_] = s2mpj_ii('LINEAR',ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(I)]
            self.A[ig,iv] = float(v_['B'+str(I)])+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('LINEAR',ig_)
        gtype = arrset(gtype,ig,'<>')
        self.gscale = arrset(self.gscale,ig,float(-0.5))
        [ig,ig_,_] = s2mpj_ii('C1',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C1')
        iv = ix_['X1']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X2']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X3']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X4']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X5']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C2',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C2')
        iv = ix_['X1']
        self.A[ig,iv] = float(10.0)+self.A[ig,iv]
        iv = ix_['X2']
        self.A[ig,iv] = float(10.0)+self.A[ig,iv]
        iv = ix_['X3']
        self.A[ig,iv] = float(-3.0)+self.A[ig,iv]
        iv = ix_['X4']
        self.A[ig,iv] = float(5.0)+self.A[ig,iv]
        iv = ix_['X5']
        self.A[ig,iv] = float(4.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C3',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C3')
        iv = ix_['X1']
        self.A[ig,iv] = float(-8.0)+self.A[ig,iv]
        iv = ix_['X2']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X3']
        self.A[ig,iv] = float(-2.0)+self.A[ig,iv]
        iv = ix_['X4']
        self.A[ig,iv] = float(-5.0)+self.A[ig,iv]
        iv = ix_['X5']
        self.A[ig,iv] = float(3.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C4',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C4')
        iv = ix_['X1']
        self.A[ig,iv] = float(8.0)+self.A[ig,iv]
        iv = ix_['X2']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X3']
        self.A[ig,iv] = float(2.0)+self.A[ig,iv]
        iv = ix_['X4']
        self.A[ig,iv] = float(5.0)+self.A[ig,iv]
        iv = ix_['X5']
        self.A[ig,iv] = float(-3.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C5',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C5')
        iv = ix_['X1']
        self.A[ig,iv] = float(-4.0)+self.A[ig,iv]
        iv = ix_['X2']
        self.A[ig,iv] = float(-2.0)+self.A[ig,iv]
        iv = ix_['X3']
        self.A[ig,iv] = float(3.0)+self.A[ig,iv]
        iv = ix_['X4']
        self.A[ig,iv] = float(-5.0)+self.A[ig,iv]
        iv = ix_['X5']
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
        self.gconst = arrset(self.gconst,ig_['C1'],float(-5.0))
        self.gconst = arrset(self.gconst,ig_['C2'],float(20.0))
        self.gconst = arrset(self.gconst,ig_['C3'],float(-40.0))
        self.gconst = arrset(self.gconst,ig_['C4'],float(11.0))
        self.gconst = arrset(self.gconst,ig_['C5'],float(-30.0))
        self.gconst = arrset(self.gconst,ig_['NONL'],float(-14463.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        self.xlower = np.zeros((self.n,1))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(1.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'ePROD', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        elftp = []
        elftp = loaset(elftp,it,0,'D')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['1']),int(v_['5'])+1):
            for J in range(int(v_['1']),int(v_['5'])+1):
                ename = 'E'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'ePROD')
                ielftype = arrset(ielftype, ie, iet_["ePROD"])
                vname = 'X'+str(I)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,1.0)
                posev = np.where(elftv[ielftype[ie]]=='X')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'X'+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,1.0)
                posev = np.where(elftv[ielftype[ie]]=='Y')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                posep = np.where(elftp[ielftype[ie]]=='D')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['D'+str(I)+','+str(J)]))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['5'])+1):
            for J in range(int(v_['1']),int(v_['5'])+1):
                ig = ig_['NONL']
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle+self.neq,self.m)] = np.zeros((self.nge,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        self.A.resize(ngrp,self.n)
        self.A     = self.A.tocsr()
        sA1,sA2    = self.A.shape
        self.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons =  np.where(self.congrps in np.setdiff1d(nlc,self.congrps))[0]
        self.pbclass = "QLR2-AN-5-5"

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def ePROD(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = self.elpar[iel_][0]*EV_[0]*EV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = self.elpar[iel_][0]*EV_[1]
            g_[1] = self.elpar[iel_][0]*EV_[0]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 0.0
                H_[1,1] = 0.0
                H_[0,1] = self.elpar[iel_][0]
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

