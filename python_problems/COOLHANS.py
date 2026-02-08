from s2mpjlib import *
class  COOLHANS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : COOLHANS
#    *********
# 
#    A problem arising from the analysis of a Cooley-Hansen economy with
#    loglinear approximation.  The problem is to solve the matrix equation
#                  A * X * X + B * X + C = 0
#    where A, B and C are known N times N matrices and X an unknown matrix
#    of matching dimension.  The instance considered here has N = 3.
# 
#    Source:
#    S. Ceria, private communication, 1995.
# 
#    SIF input: Ph. Toint, Feb 1995.
# 
#    classification = "C-CNQR2-RN-9-9"
# 
#    order of the matrix equation
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'COOLHANS'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 3
        v_['A1,1'] = 0.0
        v_['A2,1'] = 0.13725e-6
        v_['A3,1'] = 0.0
        v_['A1,2'] = 0.0
        v_['A2,2'] = 937.62
        v_['A3,2'] = 0.0
        v_['A1,3'] = 0.0
        v_['A2,3'] = -42.207
        v_['A3,3'] = 0.0
        v_['B1,1'] = 0.0060893
        v_['B2,1'] = 0.13880e-6
        v_['B3,1'] = -0.13877e-6
        v_['B1,2'] = -44.292
        v_['B2,2'] = -1886.0
        v_['B3,2'] = 42.362
        v_['B1,3'] = 2.0011
        v_['B2,3'] = 42.362
        v_['B3,3'] = -2.0705
        v_['C1,1'] = 0.0
        v_['C2,1'] = 0.0
        v_['C3,1'] = 0.0
        v_['C1,2'] = 44.792
        v_['C2,2'] = 948.21
        v_['C3,2'] = -42.684
        v_['C1,3'] = 0.0
        v_['C2,3'] = 0.0
        v_['C3,3'] = 0.0
        v_['1'] = 1
        v_['2'] = 2
        v_['3'] = 3
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['N'])+1):
            for J in range(int(v_['1']),int(v_['N'])+1):
                [iv,ix_,_] = s2mpj_ii('X'+str(I)+','+str(J),ix_)
                self.xnames=arrset(self.xnames,iv,'X'+str(I)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for K in range(int(v_['1']),int(v_['N'])+1):
            for L in range(int(v_['1']),int(v_['N'])+1):
                for M in range(int(v_['1']),int(v_['N'])+1):
                    [ig,ig_,_] = s2mpj_ii('G'+str(K)+','+str(L),ig_)
                    gtype = arrset(gtype,ig,'==')
                    cnames = arrset(cnames,ig,'G'+str(K)+','+str(L))
                    irA  = np.append(irA,[ig])
                    icA  = np.append(icA,[ix_['X'+str(M)+','+str(L)]])
                    valA = np.append(valA,float(v_['B'+str(K)+','+str(M)]))
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
        for K in range(int(v_['1']),int(v_['N'])+1):
            for L in range(int(v_['1']),int(v_['N'])+1):
                v_['-C'] = -1.0*v_['C'+str(K)+','+str(L)]
                self.gconst = arrset(self.gconst,ig_['G'+str(K)+','+str(L)],float(v_['-C']))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'XX')
        elftv = loaset(elftv,it,1,'YY')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for K in range(int(v_['1']),int(v_['N'])+1):
            for L in range(int(v_['1']),int(v_['N'])+1):
                for M in range(int(v_['1']),int(v_['N'])+1):
                    ename = 'E'+str(K)+','+str(M)+','+str(L)
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    self.elftype = arrset(self.elftype,ie,'en2PR')
                    ielftype = arrset(ielftype,ie,iet_["en2PR"])
                    self.x0 = np.zeros((self.n,1))
                    vname = 'X'+str(K)+','+str(M)
                    [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                    posev = np.where(elftv[ielftype[ie]]=='XX')[0]
                    self.elvar = loaset(self.elvar,ie,posev[0],iv)
                    vname = 'X'+str(M)+','+str(L)
                    [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                    posev = np.where(elftv[ielftype[ie]]=='YY')[0]
                    self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for L in range(int(v_['1']),int(v_['N'])+1):
            for P in range(int(v_['1']),int(v_['N'])+1):
                for M in range(int(v_['1']),int(v_['N'])+1):
                    ig = ig_['G'+str(int(v_['1']))+','+str(L)]
                    posel = len(self.grelt[ig])
                    self.grelt  = (
                          loaset(self.grelt,ig,posel,ie_['E'+str(P)+','+str(M)+','+str(L)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    self.grelw  = (
                          loaset(self.grelw,ig,posel,float(v_['A'+str(int(v_['1']))+','+str(P)])))
        for L in range(int(v_['1']),int(v_['N'])+1):
            for P in range(int(v_['1']),int(v_['N'])+1):
                for M in range(int(v_['1']),int(v_['N'])+1):
                    ig = ig_['G'+str(int(v_['2']))+','+str(L)]
                    posel = len(self.grelt[ig])
                    self.grelt  = (
                          loaset(self.grelt,ig,posel,ie_['E'+str(P)+','+str(M)+','+str(L)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    self.grelw  = (
                          loaset(self.grelw,ig,posel,float(v_['A'+str(int(v_['2']))+','+str(P)])))
        for L in range(int(v_['1']),int(v_['N'])+1):
            for P in range(int(v_['1']),int(v_['N'])+1):
                for M in range(int(v_['1']),int(v_['N'])+1):
                    ig = ig_['G'+str(int(v_['3']))+','+str(L)]
                    posel = len(self.grelt[ig])
                    self.grelt  = (
                          loaset(self.grelt,ig,posel,ie_['E'+str(P)+','+str(M)+','+str(L)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    self.grelw  = (
                          loaset(self.grelw,ig,posel,float(v_['A'+str(int(v_['3']))+','+str(P)])))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 0.0
#    Solution
# LO SOLTN                0.0
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
        self.pbclass   = "C-CNQR2-RN-9-9"
        self.x0        = np.zeros((self.n,1))
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def en2PR(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0,0]*EV_[1,0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1,0]
            g_[1] = EV_[0,0]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = 1.0
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

