from s2mpjlib import *
class  LEVYMONE6(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : LEVYMONE6
#    *********
#    A global optimization example due to Levy & Montalvo 
#    This problem is one of the parameterised set LEVYMONT5-LEVYMONT10
# 
#    Source:  problem 6 in
# 
#    A. V. Levy and A. Montalvo
#    "The Tunneling Algorithm for the Global Minimization of Functions"
#    SIAM J. Sci. Stat. Comp. 6(1) 1985 15:29 
#    https://doi.org/10.1137/0906002
# 
#    nonlinear equations version
# 
#    SIF input: Nick Gould, August 2021
# 
#    classification = "C-CNOR2-AY-3-6"
# 
#    N is the number of variables
# 
#           Alternative values for the SIF file parameters:
# IE N                   3              $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'LEVYMONE6'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(5);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
        v_['A'] = 1.0
        v_['K'] = 10.0
        v_['L'] = 0.25
        v_['C'] = 0.75
        v_['1'] = 1
        v_['2'] = 2
        v_['PI/4'] = np.arctan(1.0)
        v_['PI'] = 4.0*v_['PI/4']
        v_['RN'] = float(v_['N'])
        v_['A-C'] = v_['A']-v_['C']
        v_['PI/N'] = v_['PI']/v_['RN']
        v_['KPI/N'] = v_['K']*v_['PI/N']
        v_['ROOTKPI/N'] = np.sqrt(v_['KPI/N'])
        v_['N/PI'] = v_['RN']/v_['PI']
        v_['N/KPI'] = v_['N/PI']/v_['K']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('Q'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'Q'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(I)]])
            valA = np.append(valA,float(v_['L']))
            self.gscale = arrset(self.gscale,ig,float(v_['N/PI']))
            [ig,ig_,_] = s2mpj_ii('N'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'N'+str(I))
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
        for I in range(int(v_['1']),int(v_['N'])+1):
            self.gconst = arrset(self.gconst,ig_['Q'+str(I)],float(v_['A-C']))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-10.0)
        self.xupper = np.full((self.n,1),10.0)
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(8.0))
        self.x0[ix_['X1']] = float(-8.0)
        self.x0[ix_['X2']] = float(8.0)
        pass
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eS2', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftp = []
        elftp = loaset(elftp,it,0,'L')
        elftp = loaset(elftp,it,1,'C')
        [it,iet_,_] = s2mpj_ii( 'ePS2', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Z')
        elftp = loaset(elftp,it,0,'L')
        elftp = loaset(elftp,it,1,'C')
        elftp = loaset(elftp,it,2,'A')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        ename = 'E'+str(int(v_['1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eS2')
        ielftype = arrset(ielftype,ie,iet_["eS2"])
        ename = 'E'+str(int(v_['1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['1']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-10.0),float(10.0),float(8.0))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        posep = np.where(elftp[ielftype[ie]]=='L')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(v_['L']))
        ename = 'E'+str(int(v_['1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        posep = np.where(elftp[ielftype[ie]]=='C')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(v_['C']))
        for I in range(int(v_['2']),int(v_['N'])+1):
            v_['I-1'] = I-v_['1']
            ename = 'E'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'ePS2')
            ielftype = arrset(ielftype,ie,iet_["ePS2"])
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-10.0),float(10.0),float(8.0))
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I-1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-10.0),float(10.0),float(8.0))
            posev = np.where(elftv[ielftype[ie]]=='Z')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='L')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['L']))
            posep = np.where(elftp[ielftype[ie]]=='C')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['C']))
            posep = np.where(elftp[ielftype[ie]]=='A')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['A']))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            ig = ig_['N'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['ROOTKPI/N']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 0.0
#    Solution
# LO SOLTN               0.0
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
        self.pbclass   = "C-CNOR2-AY-3-6"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def e_globs(self):

        import numpy as np
        self.efpar = np.array([]);
        self.efpar = arrset( self.efpar,0,4.0*np.arctan(1.0e0))
        return pbm

    @staticmethod
    def eS2(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        PIL = self.efpar[0]*self.elpar[iel_][0]
        V = PIL*EV_[0,0]+self.efpar[0]*self.elpar[iel_][1]
        SINV = np.sin(V)
        COSV = np.cos(V)
        f_   = SINV
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = PIL*COSV
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = -PIL*PIL*SINV
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ePS2(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        PIL = self.efpar[0]*self.elpar[iel_][0]
        U = self.elpar[iel_][0]*EV_[1,0]+self.elpar[iel_][1]-self.elpar[iel_][2]
        V = PIL*EV_[0,0]+self.efpar[0]*self.elpar[iel_][1]
        SINV = np.sin(V)
        COSV = np.cos(V)
        f_   = U*SINV
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = PIL*U*COSV
            g_[1] = self.elpar[iel_][0]*SINV
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = -PIL*PIL*U*SINV
                H_[0,1] = self.elpar[iel_][0]*PIL*COSV
                H_[1,0] = H_[0,1]
                H_[1,1] = 0.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

