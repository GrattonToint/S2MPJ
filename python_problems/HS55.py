from s2mpjlib import *
class  HS55(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS55
#    *********
# 
#    Source: problem 55 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: A.R. Conn, April 1990
# 
#    classification = "OLR2-AN-6-6"
# 
#    some useful parameters, including N, the number of variables.
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS55'

    def __init__(self, *args): 
        import numpy as np
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 6
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
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X2']
        self.A[ig,iv] = float(2.0)+self.A[ig,iv]
        iv = ix_['X5']
        self.A[ig,iv] = float(4.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('CON1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'CON1')
        iv = ix_['X1']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X2']
        self.A[ig,iv] = float(2.0)+self.A[ig,iv]
        iv = ix_['X5']
        self.A[ig,iv] = float(5.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('CON2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'CON2')
        iv = ix_['X1']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X2']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X3']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('CON3',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'CON3')
        iv = ix_['X4']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X5']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X6']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('CON4',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'CON4')
        iv = ix_['X1']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X4']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('CON5',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'CON5')
        iv = ix_['X2']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X5']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('CON6',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'CON6')
        iv = ix_['X3']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X6']
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
        self.gconst = arrset(self.gconst,ig_['CON1'],float(6.0))
        self.gconst = arrset(self.gconst,ig_['CON2'],float(3.0))
        self.gconst = arrset(self.gconst,ig_['CON3'],float(2.0))
        self.gconst = arrset(self.gconst,ig_['CON4'],float(1.0))
        self.gconst = arrset(self.gconst,ig_['CON5'],float(2.0))
        self.gconst = arrset(self.gconst,ig_['CON6'],float(2.0))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        self.xupper[ix_['X1']] = 1.0
        self.xupper[ix_['X4']] = 1.0
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(0.0))
        if('X1' in ix_):
            self.x0[ix_['X1']] = float(1.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X1']),float(1.0)))
        if('X2' in ix_):
            self.x0[ix_['X2']] = float(2.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X2']),float(2.0)))
        if('X6' in ix_):
            self.x0[ix_['X6']] = float(2.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X6']),float(2.0)))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eEPROD', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        ename = 'E1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eEPROD')
        ielftype = arrset(ielftype, ie, iet_["eEPROD"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,0.0)
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,0.0)
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
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
        self.grelt = loaset(self.grelt,ig,posel,ie_['E1'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN                6.66666666
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
        self.pbclass = "OLR2-AN-6-6"
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eEPROD(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        PROD = np.exp(EV_[0]*EV_[1])
        f_   = PROD
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]*PROD
            g_[1] = EV_[0]*PROD
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = EV_[1]*EV_[1]*PROD
                H_[0,1] = PROD+EV_[0]*EV_[1]*PROD
                H_[1,0] = H_[0,1]
                H_[1,1] = EV_[0]*EV_[0]*PROD
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

