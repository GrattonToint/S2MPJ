from s2mpjlib import *
class  HS39(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS39
#    *********
# 
#    Source: problem 39 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: A.R. Conn, April 1990
# 
#    classification = "C-CLOR2-AN-4-2"
# 
#    some useful parameters, including N, the number of variables.
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 9 XI 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS39'

    def __init__(self, *args): 
        import numpy as np
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 4
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
        [it,iet_,_] = s2mpj_ii( 'eMSQ', iet_)
        elftv = loaset(elftv,it,0,'V1')
        [it,iet_,_] = s2mpj_ii( 'eMCB', iet_)
        elftv = loaset(elftv,it,0,'V1')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        ename = 'E1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eMCB')
        ielftype = arrset(ielftype,ie,iet_["eMCB"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(2.0))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eMSQ')
        ielftype = arrset(ielftype,ie,iet_["eMSQ"])
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(2.0))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eMSQ')
        ielftype = arrset(ielftype,ie,iet_["eMSQ"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(2.0))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eMSQ')
        ielftype = arrset(ielftype,ie,iet_["eMSQ"])
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(2.0))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
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
        self.grelt = loaset(self.grelt,ig,posel,ie_['E1'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E2'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        ig = ig_['CON2']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E3'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E4'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               -1.0
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
    def eMCB(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = -EV_[0]**3
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -3.0*EV_[0]*EV_[0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = -6.0*EV_[0]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eMSQ(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = -EV_[0]*EV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -2.0*EV_[0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = -2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

