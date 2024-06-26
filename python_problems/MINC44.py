from s2mpjlib import *
class  MINC44(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : MINC44
#    *********
# 
#    Minimize the permanent of a doubly stochastic n dimensional square matrix
#    whose trace is zero.
#    The conjecture is that the minimum is achieved when all non-diagonal
#    entries of the matrix are equal to 1/(n-1).
# 
#    Source: conjecture 44 in
#    H. Minc,
#    "Theory of Permanents 1982-1985",
#    Linear Algebra and Applications, vol. 21, pp. 109-148, 1987.
# 
#    SIF input: Ph. Toint, April 1992.
#               minor correction by Ph. Shott, Jan 1995.
# 
#    classification = "LQR2-AN-V-V"
# 
#    Size of matrix
# 
#           Alternative values for the SIF file parameters:
# IE N                   2              $-PARAMETER n = 5
# IE N                   3              $-PARAMETER n = 13
# IE N                   4              $-PARAMETER n = 27
# IE N                   5              $-PARAMETER n = 51
# IE N                   6              $-PARAMETER n = 93
# IE N                   7              $-PARAMETER n = 169
# IE N                   8              $-PARAMETER n = 311
# IE N                   9              $-PARAMETER n = 583
# IE N                   10             $-PARAMETER n = 1113
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'MINC44'

    def __init__(self, *args): 
        import numpy as np
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(5);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
        v_['0'] = 0
        v_['1'] = 1
        v_['2'] = 2
        v_['N+1'] = 1+v_['N']
        v_['N-1'] = -1+v_['N']
        v_['2**N'] = 1
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['N-I+1'] = v_['N+1']-I
            v_['I-1'] = - 1+I
            v_['R2**N'] = float(v_['2**N'])
            v_['S'+str(int(v_['N-I+1']))] = 0.1+v_['R2**N']
            v_['T'+str(int(v_['I-1']))] = 0.1+v_['R2**N']
            v_['2**N'] = v_['2**N']*v_['2']
        v_['2**N-1'] = - 1+v_['2**N']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for M in range(int(v_['1']),int(v_['N-1'])+1):
            v_['RK1'] = v_['T'+str(M)]
            v_['K1'] = int(np.fix(v_['RK1']))
            v_['K2'] = 2*v_['K1']
            v_['K1'] = 1+v_['K1']
            v_['K2'] = - 1+v_['K2']
            for K in range(int(v_['K1']),int(v_['K2'])+1):
                [iv,ix_,_] = s2mpj_ii('P'+str(K),ix_)
                self.xnames=arrset(self.xnames,iv,'P'+str(K))
        for I in range(int(v_['1']),int(v_['N'])+1):
            for J in range(int(v_['1']),int(v_['N'])+1):
                [iv,ix_,_] = s2mpj_ii('A'+str(I)+','+str(J),ix_)
                self.xnames=arrset(self.xnames,iv,'A'+str(I)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.A       = lil_matrix((1000000,1000000))
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames      = np.array([])
        self.cnames = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['P'+str(int(v_['2**N-1']))]
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        for M in range(int(v_['1']),int(v_['N-1'])+1):
            v_['RK1'] = v_['T'+str(M)]
            v_['K1'] = int(np.fix(v_['RK1']))
            v_['K2'] = 2*v_['K1']
            v_['K1'] = 1+v_['K1']
            v_['K2'] = - 1+v_['K2']
            for K in range(int(v_['K1']),int(v_['K2'])+1):
                [ig,ig_,_] = s2mpj_ii('PE'+str(K),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'PE'+str(K))
                iv = ix_['P'+str(K)]
                self.A[ig,iv] = float(- 1.0)+self.A[ig,iv]
        for I in range(int(v_['1']),int(v_['N'])+1):
            for J in range(int(v_['1']),int(v_['N'])+1):
                [ig,ig_,_] = s2mpj_ii('C'+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'C'+str(J))
                iv = ix_['A'+str(I)+','+str(J)]
                self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            for J in range(int(v_['1']),int(v_['N'])+1):
                [ig,ig_,_] = s2mpj_ii('R'+str(I),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'R'+str(I))
                iv = ix_['A'+str(I)+','+str(J)]
                self.A[ig,iv] = float(1.0)+self.A[ig,iv]
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
        for J in range(int(v_['1']),int(v_['N'])+1):
            self.gconst = arrset(self.gconst,ig_['C'+str(J)],float(1.0))
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            self.gconst = arrset(self.gconst,ig_['R'+str(I)],float(1.0))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        for I in range(int(v_['1']),int(v_['N'])+1):
            for J in range(int(v_['1']),int(v_['N'])+1):
                self.xupper[ix_['A'+str(I)+','+str(J)]] = 1.0
        for I in range(int(v_['1']),int(v_['N'])+1):
            self.xlower[ix_['A'+str(I)+','+str(I)]] = 0.0
            self.xupper[ix_['A'+str(I)+','+str(I)]] = 0.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'A')
        elftv = loaset(elftv,it,1,'P')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for M in range(int(v_['1']),int(v_['N-1'])+1):
            v_['RK1'] = v_['T'+str(M)]
            v_['K1'] = int(np.fix(v_['RK1']))
            v_['K2'] = 2*v_['K1']
            v_['K1'] = 1+v_['K1']
            v_['K2'] = - 1+v_['K2']
            for K in range(int(v_['K1']),int(v_['K2'])+1):
                v_['ID'] = 0
                v_['PT'] = 1
                v_['KK'] = K
                for I in range(int(v_['1']),int(v_['N'])+1):
                    v_['SI'] = v_['S'+str(I)]
                    v_['ISI'] = int(np.fix(v_['SI']))
                    v_['BI'] = int(np.fix(v_['KK']/v_['ISI']))
                    v_['ID'] = v_['ID']+v_['BI']
                    v_['BISI'] = v_['BI']*v_['ISI']
                    v_['KK'] = v_['KK']-v_['BISI']
                    v_['RI'] = float(I)
                    v_['RNZ'+str(int(v_['PT']))] = 0.1+v_['RI']
                    v_['PT'] = v_['PT']+v_['BI']
                v_['I1'] = v_['0']
                v_['I2'] = v_['1']
                v_['ID-2'] = - 2+v_['ID']
                for I in range(int(v_['1']),int(v_['ID-2'])+1):
                    v_['I1'] = v_['ID']
                    v_['I2'] = v_['0']
                for I in range(int(v_['1']),int(v_['I1'])+1):
                    v_['RJ'] = v_['RNZ'+str(I)]
                    v_['J'] = int(np.fix(v_['RJ']))
                    v_['SI'] = v_['S'+str(int(v_['J']))]
                    v_['ISI'] = int(np.fix(v_['SI']))
                    v_['IPP'] = K-v_['ISI']
                    ename = 'E'+str(K)+','+str(I)
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    self.elftype = arrset(self.elftype,ie,'en2PR')
                    ielftype = arrset(ielftype, ie, iet_["en2PR"])
                    vname = 'A'+str(int(v_['ID']))+','+str(int(v_['J']))
                    [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                    posev = np.where(elftv[ielftype[ie]]=='A')[0]
                    self.elvar = loaset(self.elvar,ie,posev[0],iv)
                    vname = 'P'+str(int(v_['IPP']))
                    [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                    posev = np.where(elftv[ielftype[ie]]=='P')[0]
                    self.elvar = loaset(self.elvar,ie,posev[0],iv)
                for I in range(int(v_['1']),int(v_['I2'])+1):
                    v_['RJ'] = v_['RNZ'+str(int(v_['1']))]
                    v_['J'] = int(np.fix(v_['RJ']))
                    v_['RJJ'] = v_['RNZ'+str(int(v_['2']))]
                    v_['JJ'] = int(np.fix(v_['RJJ']))
                    ename = 'E'+str(K)+','+str(int(v_['1']))
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    self.elftype = arrset(self.elftype,ie,'en2PR')
                    ielftype = arrset(ielftype, ie, iet_["en2PR"])
                    ename = 'E'+str(K)+','+str(int(v_['1']))
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    vname = 'A'+str(int(v_['2']))+','+str(int(v_['J']))
                    [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                    posev = np.where(elftv[ielftype[ie]]=='A')[0]
                    self.elvar = loaset(self.elvar,ie,posev[0],iv)
                    ename = 'E'+str(K)+','+str(int(v_['1']))
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    vname = 'A'+str(int(v_['1']))+','+str(int(v_['JJ']))
                    [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                    posev = np.where(elftv[ielftype[ie]]=='P')[0]
                    self.elvar = loaset(self.elvar,ie,posev[0],iv)
                    ename = 'E'+str(K)+','+str(int(v_['2']))
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    self.elftype = arrset(self.elftype,ie,'en2PR')
                    ielftype = arrset(ielftype, ie, iet_["en2PR"])
                    ename = 'E'+str(K)+','+str(int(v_['2']))
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    vname = 'A'+str(int(v_['2']))+','+str(int(v_['JJ']))
                    [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                    posev = np.where(elftv[ielftype[ie]]=='A')[0]
                    self.elvar = loaset(self.elvar,ie,posev[0],iv)
                    ename = 'E'+str(K)+','+str(int(v_['2']))
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    vname = 'A'+str(int(v_['1']))+','+str(int(v_['J']))
                    [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                    posev = np.where(elftv[ielftype[ie]]=='P')[0]
                    self.elvar = loaset(self.elvar,ie,posev[0],iv)
                v_['RD'] = float(v_['ID'])
                v_['D'+str(K)] = 0.1+v_['RD']
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for M in range(int(v_['1']),int(v_['N-1'])+1):
            v_['RK1'] = v_['T'+str(M)]
            v_['K1'] = int(np.fix(v_['RK1']))
            v_['K2'] = 2*v_['K1']
            v_['K1'] = 1+v_['K1']
            v_['K2'] = - 1+v_['K2']
            for K in range(int(v_['K1']),int(v_['K2'])+1):
                v_['RD'] = v_['D'+str(K)]
                v_['ID'] = int(np.fix(v_['RD']))
                for I in range(int(v_['1']),int(v_['ID'])+1):
                    ig = ig_['PE'+str(K)]
                    posel = len(self.grelt[ig])
                    self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(K)+','+str(I)])
                    nlc = np.union1d(nlc,np.array([ig]))
                    self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 0.0
#    Solution
# LO SOLTN(2)           1.0
# LO SOLTN(3)           0.25
# LO SOLTN(4)           0.11111111
# LO SOLTN(5)           0.04296835
# LO SOLTN(6)           0.01695926
# LO SOLTN(7)           6.62293832D-03
# LO SOLTN(8)           2.57309338D-03
# LO SOLTN(9)           9.94617795D-04
# LO SOLTN(10)          3.83144655D-04
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        self.A.resize(ngrp,self.n)
        self.A     = self.A.tocsr()
        sA1,sA2    = self.A.shape
        self.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons =  np.where(self.congrps in np.setdiff1d(nlc,self.congrps))[0]
        self.pbclass = "LQR2-AN-V-V"
        self.x0        = np.zeros((self.n,1))
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

