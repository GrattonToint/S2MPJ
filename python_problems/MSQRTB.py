from s2mpjlib import *
class  MSQRTB(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : MSQRTB
#    *********
# 
#    The dense matrix square root problem by Nocedal and Liu.
#    (Case 1) seen as a system of nonlinear equations.
# 
#    Source:  problem 204 (p. 93) in
#    A.R. Buckley,
#    "Test functions for unconstrained minimization",
#    TR 1989CS-3, Mathematics, statistics and computing centre,
#    Dalhousie University, Halifax (CDN), 1989.
# 
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "C-CNQR2-AN-V-V"
# 
#    Dimension of the matrix ( at least 3)
# 
#           Alternative values for the SIF file parameters:
# IE P                   3              $-PARAMETER n = 9     original value
# IE P                   7              $-PARAMETER n = 49
# IE P                   10             $-PARAMETER n = 100
# IE P                   23             $-PARAMETER n = 529
# IE P                   32             $-PARAMETER n = 1024
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'MSQRTB'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['P'] = int(5);  #  SIF file default value
        else:
            v_['P'] = int(args[0])
# IE P                   70             $-PARAMETER n = 4900
        v_['N'] = v_['P']*v_['P']
        v_['1'] = 1
        v_['K'] = 0.0
        for I in range(int(v_['1']),int(v_['P'])+1):
            for J in range(int(v_['1']),int(v_['P'])+1):
                v_['K'] = 1.0+v_['K']
                v_['K2'] = v_['K']*v_['K']
                v_['B'+str(I)+','+str(J)] = np.sin(v_['K2'])
        v_['B3,1'] = 0.0
        for I in range(int(v_['1']),int(v_['P'])+1):
            for J in range(int(v_['1']),int(v_['P'])+1):
                v_['A'+str(I)+','+str(J)] = 0.0
                for T in range(int(v_['1']),int(v_['P'])+1):
                    v_['PROD'] = v_['B'+str(I)+','+str(T)]*v_['B'+str(T)+','+str(J)]
                    v_['A'+str(I)+','+str(J)] = v_['A'+str(I)+','+str(J)]+v_['PROD']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['P'])+1):
            for J in range(int(v_['1']),int(v_['P'])+1):
                [iv,ix_,_] = s2mpj_ii('X'+str(I)+','+str(J),ix_)
                self.xnames=arrset(self.xnames,iv,'X'+str(I)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['1']),int(v_['P'])+1):
            for J in range(int(v_['1']),int(v_['P'])+1):
                [ig,ig_,_] = s2mpj_ii('G'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'G'+str(I)+','+str(J))
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
        for I in range(int(v_['1']),int(v_['P'])+1):
            for J in range(int(v_['1']),int(v_['P'])+1):
                self.gconst  = (
                      arrset(self.gconst,ig_['G'+str(I)+','+str(J)],float(v_['A'+str(I)+','+str(J)])))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        v_['K'] = 0.0
        for I in range(int(v_['1']),int(v_['P'])+1):
            for J in range(int(v_['1']),int(v_['P'])+1):
                v_['K'] = 1.0+v_['K']
                v_['K2'] = v_['K']*v_['K']
                v_['SK2'] = np.sin(v_['K2'])
                v_['-4SK2/5'] = -0.8*v_['SK2']
                v_['XIJ'] = v_['B'+str(I)+','+str(J)]+v_['-4SK2/5']
                if('X'+str(I)+','+str(J) in ix_):
                    self.x0[ix_['X'+str(I)+','+str(J)]] = float(v_['XIJ'])
                else:
                    self.y0  = (
                          arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X'+str(I)+','+str(J)]),float(v_['XIJ'])))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'XIT')
        elftv = loaset(elftv,it,1,'XTJ')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['1']),int(v_['P'])+1):
            for J in range(int(v_['1']),int(v_['P'])+1):
                for T in range(int(v_['1']),int(v_['P'])+1):
                    ename = 'E'+str(I)+','+str(J)+','+str(T)
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    self.elftype = arrset(self.elftype,ie,'en2PR')
                    ielftype = arrset(ielftype,ie,iet_["en2PR"])
                    vname = 'X'+str(I)+','+str(T)
                    [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                    posev = np.where(elftv[ielftype[ie]]=='XIT')[0]
                    self.elvar = loaset(self.elvar,ie,posev[0],iv)
                    vname = 'X'+str(T)+','+str(J)
                    [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                    posev = np.where(elftv[ielftype[ie]]=='XTJ')[0]
                    self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['P'])+1):
            for J in range(int(v_['1']),int(v_['P'])+1):
                for T in range(int(v_['1']),int(v_['P'])+1):
                    ig = ig_['G'+str(I)+','+str(J)]
                    posel = len(self.grelt[ig])
                    self.grelt  = (
                          loaset(self.grelt,ig,posel,ie_['E'+str(I)+','+str(J)+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CNQR2-AN-V-V"
        self.objderlvl = 2
        self.conderlvl = [2]


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

