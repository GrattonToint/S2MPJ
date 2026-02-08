from s2mpjlib import *
class  HONG(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    Source: Se June Hong/Chid Apte
# 
#    SIF input: A.R.Conn, Jan 1991.
# 
#    classification = "C-COLR2-AN-4-1"
# 
#   Problem parameters
# 
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HONG'

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
        [iv,ix_,_] = s2mpj_ii('T1',ix_)
        self.xnames=arrset(self.xnames,iv,'T1')
        [iv,ix_,_] = s2mpj_ii('T2',ix_)
        self.xnames=arrset(self.xnames,iv,'T2')
        [iv,ix_,_] = s2mpj_ii('T3',ix_)
        self.xnames=arrset(self.xnames,iv,'T3')
        [iv,ix_,_] = s2mpj_ii('T4',ix_)
        self.xnames=arrset(self.xnames,iv,'T4')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('SUM1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'SUM1')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['T1']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['T2']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['T3']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['T4']])
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
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        self.gconst = arrset(self.gconst,ig_['SUM1'],float(1.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),0.0)
        self.xupper = np.full((self.n,1),1.0)
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(0.5))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eEXP', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftp = []
        elftp = loaset(elftp,it,0,'P1')
        elftp = loaset(elftp,it,1,'P2')
        elftp = loaset(elftp,it,2,'P3')
        elftp = loaset(elftp,it,3,'P4')
        elftp = loaset(elftp,it,4,'P5')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        ename = 'E1'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eEXP')
            ielftype = arrset(ielftype,ie,iet_['eEXP'])
        vname = 'T1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.0),float(1.0),float(0.5))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='P1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(0.0))
        posep = np.where(elftp[ielftype[ie]]=='P2')[0]
        loaset(self.elpar,ie,posep[0],float(25.0))
        posep = np.where(elftp[ielftype[ie]]=='P3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(0.92))
        posep = np.where(elftp[ielftype[ie]]=='P4')[0]
        loaset(self.elpar,ie,posep[0],float(0.08))
        posep = np.where(elftp[ielftype[ie]]=='P5')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(0.38))
        ename = 'E2'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eEXP')
            ielftype = arrset(ielftype,ie,iet_['eEXP'])
        vname = 'T2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.0),float(1.0),float(0.5))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='P1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(0.0))
        posep = np.where(elftp[ielftype[ie]]=='P2')[0]
        loaset(self.elpar,ie,posep[0],float(50.0))
        posep = np.where(elftp[ielftype[ie]]=='P3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-2.95))
        posep = np.where(elftp[ielftype[ie]]=='P4')[0]
        loaset(self.elpar,ie,posep[0],float(3.95))
        posep = np.where(elftp[ielftype[ie]]=='P5')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(0.11))
        ename = 'E3'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eEXP')
            ielftype = arrset(ielftype,ie,iet_['eEXP'])
        vname = 'T3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.0),float(1.0),float(0.5))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='P1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(9.0))
        posep = np.where(elftp[ielftype[ie]]=='P2')[0]
        loaset(self.elpar,ie,posep[0],float(-4.0))
        posep = np.where(elftp[ielftype[ie]]=='P3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.66))
        posep = np.where(elftp[ielftype[ie]]=='P4')[0]
        loaset(self.elpar,ie,posep[0],float(1657834.0))
        posep = np.where(elftp[ielftype[ie]]=='P5')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.48))
        ename = 'E4'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eEXP')
            ielftype = arrset(ielftype,ie,iet_['eEXP'])
        vname = 'T4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.0),float(1.0),float(0.5))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='P1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(0.0))
        posep = np.where(elftp[ielftype[ie]]=='P2')[0]
        loaset(self.elpar,ie,posep[0],float(20000.0))
        posep = np.where(elftp[ielftype[ie]]=='P3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(0.11))
        posep = np.where(elftp[ielftype[ie]]=='P4')[0]
        loaset(self.elpar,ie,posep[0],float(0.89))
        posep = np.where(elftp[ielftype[ie]]=='P5')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(0.00035))
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
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E2'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E3'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E4'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution unknown
        self.objlower = -4.0
        self.objupper = 300.0
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-COLR2-AN-4-1"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eEXP(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        XTOT = self.elpar[iel_][0]+self.elpar[iel_][1]*EV_[0,0]
        EP5 = np.exp(self.elpar[iel_][4]*XTOT)
        f_   = self.elpar[iel_][2]+self.elpar[iel_][3]*EP5
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = self.elpar[iel_][1]*self.elpar[iel_][3]*self.elpar[iel_][4]*EP5
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0]  = (
                      self.elpar[iel_][1]*self.elpar[iel_][1]*self.elpar[iel_][3]*self.elpar[iel_][4]*self.elpar[iel_][4]*EP5)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

