from s2mpjlib import *
class  ELATTAR(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : ELATTAR
#    *********
# 
#    A nonlinear minmax problem in six variables.
# 
#    The problem is nonconvex and has several local minima.
# 
#    Source: 
#    R.A. El-Attar, M. Vidyasagar and S.R.K. Dutta,
#    "An algorithm for l_1-approximation",
#    SINUM 16, pp.70-86, 1979.
# 
#    SIF input: Ph. Toint, Nov 1993.
# 
#    classification = "C-CLOR2-AN-7-102"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'ELATTAR'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['1'] = 1
        v_['6'] = 6
        v_['51'] = 51
        v_['T'] = 0.0
        for I in range(int(v_['1']),int(v_['51'])+1):
            v_['T'+str(I)] = v_['T']
            v_['T'] = 0.1+v_['T']
            v_['ETI'] = np.exp(v_['T'+str(I)])
            v_['Y'+str(I)] = 0.5*v_['ETI']
            v_['-2TI'] = -2.0*v_['T'+str(I)]
            v_['E-2TI'] = np.exp(v_['-2TI'])
            v_['Y'+str(I)] = v_['Y'+str(I)]-v_['E-2TI']
            v_['-3TI'] = -3.0*v_['T'+str(I)]
            v_['E-3TI'] = np.exp(v_['-3TI'])
            v_['E-3TI/2'] = 0.5*v_['E-3TI']
            v_['Y'+str(I)] = v_['Y'+str(I)]+v_['E-3TI/2']
            v_['-3TI/2'] = 0.5*v_['-3TI']
            v_['E-3TI/2'] = np.exp(v_['-3TI/2'])
            v_['7TI'] = 7.0*v_['T'+str(I)]
            v_['S7TI'] = np.sin(v_['7TI'])
            v_['TT'] = v_['E-3TI/2']*v_['S7TI']
            v_['TT'] = 1.5*v_['TT']
            v_['Y'+str(I)] = v_['Y'+str(I)]+v_['TT']
            v_['5TI'] = 5.0*v_['T'+str(I)]
            v_['-5TI/2'] = -0.5*v_['5TI']
            v_['E-5TI/2'] = np.exp(v_['-5TI/2'])
            v_['S5TI'] = np.sin(v_['5TI'])
            v_['TT'] = v_['E-5TI/2']*v_['S5TI']
            v_['Y'+str(I)] = v_['Y'+str(I)]+v_['TT']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['6'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I))
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
        for I in range(int(v_['1']),int(v_['51'])+1):
            [ig,ig_,_] = s2mpj_ii('F'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'F'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['U']])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('MF'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'MF'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['U']])
            valA = np.append(valA,float(-1.0))
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
        for I in range(int(v_['1']),int(v_['51'])+1):
            self.gconst = arrset(self.gconst,ig_['F'+str(I)],float(v_['Y'+str(I)]))
            v_['-Y'+str(I)] = -1.0*v_['Y'+str(I)]
            self.gconst = arrset(self.gconst,ig_['MF'+str(I)],float(v_['-Y'+str(I)]))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        if('X1' in ix_):
            self.x0[ix_['X1']] = float(-2.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X1']),float(-2.0)))
        if('X2' in ix_):
            self.x0[ix_['X2']] = float(-2.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X2']),float(-2.0)))
        if('X3' in ix_):
            self.x0[ix_['X3']] = float(7.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X3']),float(7.0)))
        if('X5' in ix_):
            self.x0[ix_['X5']] = float(-2.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X5']),float(-2.0)))
        if('X6' in ix_):
            self.x0[ix_['X6']] = float(1.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X6']),float(1.0)))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eET1', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftv = loaset(elftv,it,3,'V4')
        elftp = []
        elftp = loaset(elftp,it,0,'T')
        [it,iet_,_] = s2mpj_ii( 'eET2', iet_)
        elftv = loaset(elftv,it,0,'V5')
        elftv = loaset(elftv,it,1,'V6')
        elftp = loaset(elftp,it,0,'T')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['1']),int(v_['51'])+1):
            ename = 'EL1'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eET1')
            ielftype = arrset(ielftype,ie,iet_["eET1"])
            vname = 'X1'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X4'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V4')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='T')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['T'+str(I)]))
            ename = 'EL2'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eET2')
            ielftype = arrset(ielftype,ie,iet_["eET2"])
            vname = 'X5'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V5')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X6'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V6')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='T')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['T'+str(I)]))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['51'])+1):
            ig = ig_['F'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['EL1'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(1.0))
            posel = posel+1
            self.grelt = loaset(self.grelt,ig,posel,ie_['EL2'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,float(1.0))
            ig = ig_['MF'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['EL1'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
            posel = posel+1
            self.grelt = loaset(self.grelt,ig,posel,ie_['EL2'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 0.0
#    Solution       
# LO SOLTN               0.1427066255
# LO SOLTN               74.206179244
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.cupper[np.arange(self.nle)] = np.zeros((self.nle,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CLOR2-AN-7-102"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eET1(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        A = -EV_[1]*self.elpar[iel_][0]
        B = EV_[2]*self.elpar[iel_][0]+EV_[3]
        EA = np.exp(A)
        CB = np.cos(B)
        SB = np.sin(B)
        EACB = EA*CB
        EASB = EA*SB
        V1EACB = EV_[0]*EACB
        V1EASB = EV_[0]*EASB
        f_   = V1EACB
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EACB
            g_[1] = -self.elpar[iel_][0]*V1EACB
            g_[2] = -self.elpar[iel_][0]*V1EASB
            g_[3] = -V1EASB
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,1] = -self.elpar[iel_][0]*EACB
                H_[1,0] = H_[0,1]
                H_[0,2] = -self.elpar[iel_][0]*EASB
                H_[2,0] = H_[0,2]
                H_[0,3] = -EASB
                H_[3,0] = H_[0,3]
                H_[1,1] = self.elpar[iel_][0]*self.elpar[iel_][0]*V1EACB
                H_[1,2] = self.elpar[iel_][0]*self.elpar[iel_][0]*V1EASB
                H_[2,1] = H_[1,2]
                H_[1,3] = self.elpar[iel_][0]*V1EASB
                H_[3,1] = H_[1,3]
                H_[2,2] = -self.elpar[iel_][0]*self.elpar[iel_][0]*V1EACB
                H_[2,3] = -self.elpar[iel_][0]*V1EACB
                H_[3,2] = H_[2,3]
                H_[3,3] = -V1EACB
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eET2(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        A = -EV_[1]*self.elpar[iel_][0]
        EA = np.exp(A)
        B = EV_[0]*EA
        f_   = B
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EA
            g_[1] = -self.elpar[iel_][0]*B
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = -self.elpar[iel_][0]*EA
                H_[1,0] = H_[0,1]
                H_[1,1] = self.elpar[iel_][0]*self.elpar[iel_][0]*B
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

