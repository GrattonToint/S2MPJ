from s2mpjlib import *
class  HS64(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS64
#    *********
# 
#    Source: problem 64 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: N. Gould, Dec 1989.
# 
#    classification = "C-COOR2-AN-3-1"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 9 XI 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS64'

    def __init__(self, *args): 
        import numpy as np
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['1'] = 1
        v_['3'] = 3
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
        iv = ix_['X1']
        self.A[ig,iv] = float(5.0)+self.A[ig,iv]
        iv = ix_['X2']
        self.A[ig,iv] = float(20.0)+self.A[ig,iv]
        iv = ix_['X3']
        self.A[ig,iv] = float(10.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('CONSTR',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'CONSTR')
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
        self.gconst = arrset(self.gconst,ig_['CONSTR'],float(1.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),0.00001)
        self.xupper = np.full((self.n,1),+float('inf'))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(1.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eRECIP', iet_)
        elftv = loaset(elftv,it,0,'V1')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['1']),int(v_['3'])+1):
            ename = 'OBJE'+str(I)
            [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
            if newelt:
                self.elftype = arrset(self.elftype,ie,'eRECIP')
                ielftype = arrset(ielftype,ie,iet_['eRECIP'])
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.00001),None,float(1.0))
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'CONE'+str(I)
            [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eRECIP')
            ielftype = arrset(ielftype,ie,iet_["eRECIP"])
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.00001),None,float(1.0))
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
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
        self.grelt = loaset(self.grelt,ig,posel,ie_['OBJE1'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(50000.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['OBJE2'])
        self.grelw = loaset(self.grelw,ig,posel,float(72000.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['OBJE3'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(144000.0))
        ig = ig_['CONSTR']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['CONE1'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(4.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['CONE2'])
        self.grelw = loaset(self.grelw,ig,posel,float(32.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['CONE3'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(120.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               6299.842428
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
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
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass = "C-COOR2-AN-3-1"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eRECIP(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = 1.0/EV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -1.0/(EV_[0]**2)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0/(EV_[0]**3)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

