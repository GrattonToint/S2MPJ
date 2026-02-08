from s2mpjlib import *
class  TWOBARS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    Structureal analysis of the simplest two bar scheme.  The structure has
#    the following simple symmetric shape
# 
#                                 *
#                                / \
#                               /   \
#                              /     \
#                            """     """
# 
#    and a force is applied at the top node.  The unknown are the distance
#    of the left and right feet wrt to the projection of the top node and the
#    weight of the bars.
# 
#    Source:
#    an example in a talk by W.K. Zhang and C. Fleury, LLN, 1994.
# 
#    SIF input: Ph. Toint, November 1994
# 
#    classification = "C-COOR2-MN-2-2"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'TWOBARS'

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
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('CONS1',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'CONS1')
        [ig,ig_,_] = s2mpj_ii('CONS2',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'CONS2')
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
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        self.gconst = arrset(self.gconst,ig_['CONS1'],float(1.0))
        self.gconst = arrset(self.gconst,ig_['CONS2'],float(1.0))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        self.xlower[ix_['X1']] = 0.2
        self.xupper[ix_['X1']] = 4.0
        self.xlower[ix_['X2']] = 0.1
        self.xupper[ix_['X2']] = 1.6
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(1.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eOE', iet_)
        elftv = loaset(elftv,it,0,'XX')
        elftv = loaset(elftv,it,1,'YY')
        [it,iet_,_] = s2mpj_ii( 'eCE1', iet_)
        elftv = loaset(elftv,it,0,'XX')
        elftv = loaset(elftv,it,1,'YY')
        [it,iet_,_] = s2mpj_ii( 'eCE2', iet_)
        elftv = loaset(elftv,it,0,'XX')
        elftv = loaset(elftv,it,1,'YY')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        ename = 'OBEL'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eOE')
        ielftype = arrset(ielftype,ie,iet_["eOE"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(1.0))
        posev = np.where(elftv[ielftype[ie]]=='XX')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(1.0))
        posev = np.where(elftv[ielftype[ie]]=='YY')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'COEL1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCE1')
        ielftype = arrset(ielftype,ie,iet_["eCE1"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(1.0))
        posev = np.where(elftv[ielftype[ie]]=='XX')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(1.0))
        posev = np.where(elftv[ielftype[ie]]=='YY')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'COEL2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCE2')
        ielftype = arrset(ielftype,ie,iet_["eCE2"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(1.0))
        posev = np.where(elftv[ielftype[ie]]=='XX')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(1.0))
        posev = np.where(elftv[ielftype[ie]]=='YY')[0]
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
        self.grelt = loaset(self.grelt,ig,posel,ie_['OBEL'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['CONS1']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['COEL1'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.124))
        ig = ig_['CONS2']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['COEL2'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.124))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 0.0
#    Solution
# LO SOLTN               1.5086379655
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.cupper[np.arange(self.nle)] = np.zeros((self.nle,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-COOR2-MN-2-2"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eOE(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        A = 1.0+EV_[1,0]*EV_[1,0]
        RA = np.sqrt(A)
        f_   = EV_[0,0]*RA
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = RA
            g_[1] = EV_[0,0]*EV_[1,0]/RA
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = EV_[1,0]/RA
                H_[1,0] = H_[0,1]
                H_[1,1] = EV_[0,0]/(A*RA)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eCE1(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        A = 1.0+EV_[1,0]*EV_[1,0]
        RA = np.sqrt(A)
        B = 8.0/EV_[0,0]
        DB = -8.0/EV_[0,0]**2
        D2B = 16.0/EV_[0,0]**3
        C = 1.0/(EV_[0,0]*EV_[1,0])
        DCDX = -1.0/(EV_[0,0]**2*EV_[1,0])
        DCDY = -1.0/(EV_[1,0]**2*EV_[0,0])
        D2CDXX = 2.0/(EV_[0,0]**3*EV_[1,0])
        D2CDXY = 1.0/(EV_[0,0]*EV_[1,0])**2
        D2CDYY = 2.0/(EV_[0,0]*EV_[1,0]**3)
        BC = B+C
        f_   = RA*BC
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = RA*(DB+DCDX)
            g_[1] = EV_[1,0]*BC/RA+RA*DCDY
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = RA*(D2B+D2CDXX)
                H_[0,1] = RA*D2CDXY+EV_[1,0]*(DB+DCDX)/RA
                H_[1,0] = H_[0,1]
                H_[1,1] = (BC+2.0*EV_[1,0]*DCDY-EV_[1,0]*EV_[1,0]*BC/A)/RA+RA*D2CDYY
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eCE2(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        A = 1.0+EV_[1,0]*EV_[1,0]
        RA = np.sqrt(A)
        B = 8.0/EV_[0,0]
        DB = -8.0/EV_[0,0]**2
        D2B = 16.0/EV_[0,0]**3
        C = 1.0/(EV_[0,0]*EV_[1,0])
        DCDX = -1.0/(EV_[0,0]**2*EV_[1,0])
        DCDY = -1.0/(EV_[1,0]**2*EV_[0,0])
        D2CDXX = 2.0/(EV_[0,0]**3*EV_[1,0])
        D2CDXY = 1.0/(EV_[0,0]*EV_[1,0])**2
        D2CDYY = 2.0/(EV_[0,0]*EV_[1,0]**3)
        BC = B-C
        f_   = RA*BC
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = RA*(DB-DCDX)
            g_[1] = EV_[1,0]*BC/RA-RA*DCDY
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = RA*(D2B-D2CDXX)
                H_[0,1] = -RA*D2CDXY+EV_[1,0]*(DB-DCDX)/RA
                H_[1,0] = H_[0,1]
                H_[1,1] = (BC-2.0*EV_[1,0]*DCDY-EV_[1,0]*EV_[1,0]*BC/A)/RA-RA*D2CDYY
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

