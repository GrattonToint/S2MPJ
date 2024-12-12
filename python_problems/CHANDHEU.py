from s2mpjlib import *
class  CHANDHEU(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CHANDHEU
#    *********
#    Chandrasekhar Radiative Transfer H equation, as stated by T. Kelley.
# 
#    This is an unconstrained version of problem CHANDHEQ.
# 
#    Source: problem 4 in
#    J.J. More',
#    "A collection of nonlinear model problems"
#    Proceedings of the AMS-SIAM Summer seminar on the Computational
#    Solution of Nonlinear Systems of Equations, Colorado, 1988.
#    Argonne National Laboratory MCS-P60-0289, 1989.
# 
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "C-CNOR2-RN-V-V"
# 
#    N is the number of discretization points
# 
#           Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER     original value
# IE N                   50             $-PARAMETER
# IE N                   100            $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 25 XI 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'CHANDHEU'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(500);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
        if nargin<2:
            v_['C'] = float(1.0);  #  SIF file default value
        else:
            v_['C'] = float(args[1])
        v_['1'] = 1
        v_['1.0'] = 1.0
        v_['RN'] = float(v_['N'])
        v_['1/N'] = v_['1.0']/v_['RN']
        v_['HALFC'] = 0.5*v_['C']
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['RI'] = float(I)
            v_['X'+str(I)] = v_['RI']*v_['1/N']
            v_['W'+str(I)] = 0.0+v_['1/N']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('H'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'H'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('G'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'G'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['H'+str(I)]])
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
        #%%%%%%%%%%%%%%%%%%  CONSTANTS %%%%%%%%%%%%%%%%%%%
        self.gconst = np.full((ngrp,1),1.0)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(1.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'HI')
        elftv = loaset(elftv,it,1,'HJ')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['1']),int(v_['N'])+1):
            for J in range(int(v_['1']),int(v_['N'])+1):
                ename = 'A'+str(I)+','+str(J)
                [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
                if newelt:
                    self.elftype = arrset(self.elftype,ie,'en2PR')
                    ielftype = arrset(ielftype,ie,iet_['en2PR'])
                vname = 'H'+str(I)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(1.0))
                posev = np.where(elftv[ielftype[ie]]=='HI')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'H'+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(1.0))
                posev = np.where(elftv[ielftype[ie]]=='HJ')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['T1'] = v_['HALFC']*v_['X'+str(I)]
            for J in range(int(v_['1']),int(v_['N'])+1):
                v_['T3'] = v_['T1']*v_['W'+str(J)]
                v_['XIPXJ'] = v_['X'+str(I)]+v_['X'+str(J)]
                v_['VAL'] = v_['T3']/v_['XIPXJ']
                v_['COEFF'] = -1.0*v_['VAL']
                ig = ig_['G'+str(I)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['A'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(v_['COEFF']))
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
        self.pbclass   = "C-CNOR2-RN-V-V"
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
        f_   = EV_[0]*EV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]
            g_[1] = EV_[0]
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

