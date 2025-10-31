from s2mpjlib import *
class  STNQP2(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : STNQP2
#    *********
# 
#    Another non-convex quadratic program with some structure.
# 
#    The objective function is of the form
#       sum (i=0,n) x_i^2 - 0.5 sum (l=1,n/p) sum(i=1,p) sum(k;i) x_{k+l}^2,
#    where n = 2^p and (k;i) means k takes the values of the first i powers of 2
#    eg, (k:3) = {k = {1,2,4}} and (k:7) = {k = {1,2,4,8,16,32}}.
#    There are equality constraints of the form
#    
#       sum(j=1,i) x_{(l-1)p+i} = i, where l=1,n/p,2 and i=1,p.
#    Finally, there are simple bounds
#          2 <= x_i, y_i <= 2    (i=0,n).
# 
#    SIF input: Nick Gould, May 1996
# 
#    classification = "C-CQLR2-AN-V-V"
# 
#    There will be 2**p + 1 variables
# 
#           Alternative values for the SIF file parameters:
# IE P                   2              $-PARAMETER n = 5
# IE P                   4              $-PARAMETER n = 17
# IE P                   6              $-PARAMETER n = 65
# IE P                   8              $-PARAMETER n = 257
# IE P                   10             $-PARAMETER n = 1025
# IE P                   12             $-PARAMETER n = 4097     original value
# IE P                   13             $-PARAMETER n = 8193
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'STNQP2'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['P'] = int(4);  #  SIF file default value
        else:
            v_['P'] = int(args[0])
# IE P                   14             $-PARAMETER n = 16395
# IE P                   15             $-PARAMETER n = 32769
# IE P                   16             $-PARAMETER n = 65537
        v_['0'] = 0
        v_['1'] = 1
        v_['2'] = 2
        v_['N'] = 1
        for I in range(int(v_['1']),int(v_['P'])+1):
            v_['N'] = v_['N']*v_['2']
        v_['N/P'] = int(np.fix(v_['N']/v_['P']))
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['0']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['0']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('O'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(I)]])
            valA = np.append(valA,float(1.0))
        for L in range(int(v_['1']),int(v_['N/P'])+1):
            for I in range(int(v_['1']),int(v_['P'])+1):
                v_['K'] = v_['1']
                for J in range(int(v_['1']),int(I)+1):
                    v_['K+L'] = v_['K']+L
                    [ig,ig_,_] = s2mpj_ii('N'+str(I)+','+str(L),ig_)
                    gtype = arrset(gtype,ig,'<>')
                    irA  = np.append(irA,[ig])
                    icA  = np.append(icA,[ix_['X'+str(int(v_['K+L']))]])
                    valA = np.append(valA,float(1.0))
                    v_['K'] = v_['K']*v_['2']
        for L in range(int(v_['1']),int(v_['N/P'])+1,int(v_['2'])):
            v_['LL'] = L*v_['P']
            v_['LL'] = v_['LL']-v_['P']
            for I in range(int(v_['1']),int(v_['P'])+1):
                for J in range(int(v_['1']),int(I)+1):
                    v_['LL+J'] = v_['LL']+J
                    [ig,ig_,_] = s2mpj_ii('E'+str(I)+','+str(L),ig_)
                    gtype = arrset(gtype,ig,'==')
                    cnames = arrset(cnames,ig,'E'+str(I)+','+str(L))
                    irA  = np.append(irA,[ig])
                    icA  = np.append(icA,[ix_['X'+str(int(v_['LL+J']))]])
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
        for L in range(int(v_['1']),int(v_['N/P'])+1,int(v_['2'])):
            for I in range(int(v_['1']),int(v_['P'])+1):
                v_['RI'] = float(I)
                self.gconst = arrset(self.gconst,ig_['E'+str(I)+','+str(L)],float(v_['RI']))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        for I in range(int(v_['0']),int(v_['N'])+1):
            self.xlower[ix_['X'+str(I)]] = -2.0
            self.xupper[ix_['X'+str(I)]] = 2.0
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(0.5))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gPSQR',igt_)
        [it,igt_,_] = s2mpj_ii('gPSQR',igt_)
        grftp = []
        grftp = loaset(grftp,it,0,'P')
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        self.grpar   = []
        for I in range(int(v_['0']),int(v_['N'])+1):
            ig = ig_['O'+str(I)]
            self.grftype = arrset(self.grftype,ig,'gPSQR')
            posgp = np.where(grftp[igt_[self.grftype[ig]]]=='P')[0]
            self.grpar =loaset(self.grpar,ig,posgp[0],float(1.0))
        for L in range(int(v_['1']),int(v_['N/P'])+1):
            for I in range(int(v_['1']),int(v_['P'])+1):
                ig = ig_['N'+str(I)+','+str(L)]
                self.grftype = arrset(self.grftype,ig,'gPSQR')
                posgp = np.where(grftp[igt_[self.grftype[ig]]]=='P')[0]
                self.grpar =loaset(self.grpar,ig,posgp[0],float(-0.5))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLUTION            -2.476395E+5   $ (P=12)
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons   = np.arange(len(self.congrps))
        self.pbclass   = "C-CQLR2-AN-V-V"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gPSQR(self,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= self.grpar[igr_][0]*GVAR_*GVAR_
        if nargout>1:
            g_ = 2.0*self.grpar[igr_][0]*GVAR_
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 2.0*self.grpar[igr_][0]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

