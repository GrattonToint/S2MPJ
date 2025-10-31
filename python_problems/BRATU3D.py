from s2mpjlib import *
class  BRATU3D(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : BRATU3D
#    *********
# 
#    The 3D Bratu problem on the unit cube, using finite differences,
# 
#    Source: Problem 3 in
#    J.J. More',
#    "A collection of nonlinear model problems"
#    Proceedings of the AMS-SIAM Summer seminar on the Computational
#    Solution of Nonlinear Systems of Equations, Colorado, 1988.
#    Argonne National Laboratory MCS-P60-0289, 1989.
# 
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "C-CNOR2-MN-V-V"
# 
#    P is the number of points in one side of the unit cube
#    The number of variables is equal to P**3
# 
#           Alternative values for the SIF file parameters:
# IE P                   3              $-PARAMETER  n = 27   original value
# IE P                   5              $-PARAMETER  n = 125
# IE P                   8              $-PARAMETER  n = 512
# IE P                   10             $-PARAMETER  n = 1000
# IE P                   17             $-PARAMETER  n = 4913
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'BRATU3D'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['P'] = int(3);  #  SIF file default value
        else:
            v_['P'] = int(args[0])
        if nargin<2:
            v_['LAMBDA'] = float(6.80812);  #  SIF file default value
        else:
            v_['LAMBDA'] = float(args[1])
        v_['1'] = 1
        v_['2'] = 2
        v_['1.0'] = 1.0
        v_['P-1'] = -1+v_['P']
        v_['RP-1'] = float(v_['P-1'])
        v_['H'] = v_['1.0']/v_['RP-1']
        v_['H2'] = v_['H']*v_['H']
        v_['C'] = v_['H2']*v_['LAMBDA']
        v_['-C'] = -1.0*v_['C']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for J in range(int(v_['1']),int(v_['P'])+1):
            for I in range(int(v_['1']),int(v_['P'])+1):
                for K in range(int(v_['1']),int(v_['P'])+1):
                    [iv,ix_,_] = s2mpj_ii('U'+str(I)+','+str(J)+','+str(K),ix_)
                    self.xnames=arrset(self.xnames,iv,'U'+str(I)+','+str(J)+','+str(K))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['2']),int(v_['P-1'])+1):
            v_['R'] = 1+I
            v_['S'] = -1+I
            for J in range(int(v_['2']),int(v_['P-1'])+1):
                v_['V'] = 1+J
                v_['W'] = -1+J
                for K in range(int(v_['2']),int(v_['P-1'])+1):
                    v_['Y'] = 1+K
                    v_['Z'] = -1+K
                    [ig,ig_,_] = s2mpj_ii('G'+str(I)+','+str(J)+','+str(K),ig_)
                    gtype = arrset(gtype,ig,'==')
                    cnames = arrset(cnames,ig,'G'+str(I)+','+str(J)+','+str(K))
                    irA  = np.append(irA,[ig])
                    icA  = np.append(icA,[ix_['U'+str(I)+','+str(J)+','+str(K)]])
                    valA = np.append(valA,float(6.0))
                    irA  = np.append(irA,[ig])
                    icA  = np.append(icA,[ix_['U'+str(int(v_['R']))+','+str(J)+','+str(K)]])
                    valA = np.append(valA,float(-1.0))
                    irA  = np.append(irA,[ig])
                    icA  = np.append(icA,[ix_['U'+str(int(v_['S']))+','+str(J)+','+str(K)]])
                    valA = np.append(valA,float(-1.0))
                    irA  = np.append(irA,[ig])
                    icA  = np.append(icA,[ix_['U'+str(I)+','+str(int(v_['V']))+','+str(K)]])
                    valA = np.append(valA,float(-1.0))
                    irA  = np.append(irA,[ig])
                    icA  = np.append(icA,[ix_['U'+str(I)+','+str(int(v_['W']))+','+str(K)]])
                    valA = np.append(valA,float(-1.0))
                    irA  = np.append(irA,[ig])
                    icA  = np.append(icA,[ix_['U'+str(I)+','+str(J)+','+str(int(v_['Y']))]])
                    valA = np.append(valA,float(-1.0))
                    irA  = np.append(irA,[ig])
                    icA  = np.append(icA,[ix_['U'+str(I)+','+str(J)+','+str(int(v_['Z']))]])
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
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        for J in range(int(v_['1']),int(v_['P'])+1):
            for K in range(int(v_['1']),int(v_['P'])+1):
                self.xlower[ix_['U'+str(int(v_['1']))+','+str(J)+','+str(K)]] = 0.0
                self.xupper[ix_['U'+str(int(v_['1']))+','+str(J)+','+str(K)]] = 0.0
                self.xlower[ix_['U'+str(int(v_['P']))+','+str(J)+','+str(K)]] = 0.0
                self.xupper[ix_['U'+str(int(v_['P']))+','+str(J)+','+str(K)]] = 0.0
        for I in range(int(v_['2']),int(v_['P-1'])+1):
            for K in range(int(v_['1']),int(v_['P'])+1):
                self.xlower[ix_['U'+str(I)+','+str(int(v_['P']))+','+str(K)]] = 0.0
                self.xupper[ix_['U'+str(I)+','+str(int(v_['P']))+','+str(K)]] = 0.0
                self.xlower[ix_['U'+str(I)+','+str(int(v_['1']))+','+str(K)]] = 0.0
                self.xupper[ix_['U'+str(I)+','+str(int(v_['1']))+','+str(K)]] = 0.0
        for I in range(int(v_['2']),int(v_['P-1'])+1):
            for J in range(int(v_['2']),int(v_['P-1'])+1):
                self.xlower[ix_['U'+str(I)+','+str(J)+','+str(int(v_['1']))]] = 0.0
                self.xupper[ix_['U'+str(I)+','+str(J)+','+str(int(v_['1']))]] = 0.0
                self.xlower[ix_['U'+str(I)+','+str(J)+','+str(int(v_['P']))]] = 0.0
                self.xupper[ix_['U'+str(I)+','+str(J)+','+str(int(v_['P']))]] = 0.0
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(0.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eEXP', iet_)
        elftv = loaset(elftv,it,0,'U')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['2']),int(v_['P-1'])+1):
            for J in range(int(v_['2']),int(v_['P-1'])+1):
                for K in range(int(v_['2']),int(v_['P-1'])+1):
                    ename = 'A'+str(I)+','+str(J)+','+str(K)
                    [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
                    if newelt:
                        self.elftype = arrset(self.elftype,ie,'eEXP')
                        ielftype = arrset(ielftype,ie,iet_['eEXP'])
                    vname = 'U'+str(I)+','+str(J)+','+str(K)
                    [iv,ix_]  = (
                          s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0)))
                    posev = np.where(elftv[ielftype[ie]]=='U')[0]
                    self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['2']),int(v_['P-1'])+1):
            for J in range(int(v_['2']),int(v_['P-1'])+1):
                for K in range(int(v_['2']),int(v_['P-1'])+1):
                    ig = ig_['G'+str(I)+','+str(J)+','+str(K)]
                    posel = len(self.grelt[ig])
                    self.grelt  = (
                          loaset(self.grelt,ig,posel,ie_['A'+str(I)+','+str(J)+','+str(K)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    self.grelw = loaset(self.grelw,ig,posel,float(v_['-C']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 0.0
#    Solution
        pass
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
        self.pbclass   = "C-CNOR2-MN-V-V"
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
        EXPU = np.exp(EV_[0])
        f_   = EXPU
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EXPU
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = EXPU
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

