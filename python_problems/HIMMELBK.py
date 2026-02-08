from s2mpjlib import *
class  HIMMELBK(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HIMMELBK
#    *********
# 
#    A blending problem for multi-component mixtures, by Paviani.
#    It has a linear objective and linear and nonlinear constraints.
# 
#    Compared to the problem specified in Himmelblau, the inequality
#    constraints have been removed, because, as stated in this source,
#    they impose that
#    X(1)=X(2)=X(3)=X(7)=X(9)=X(9)=X(13)=X(14)=X(15)=X(19)=X(20)=X(21)=0
#    which is clearly contradictory with the given solution(s).  As
#    there does not seem to be a natural way to correct this statement
#    without knowing more about the original problem, the troublesome
#    constraints have been removed.
# 
#    Source: from problem 20 in
#    D.H. Himmelblau,
#    "Applied nonlinear programming",
#    McGraw-Hill, New-York, 1972.
# 
#    SIF input: Ph. Toint, March 1991.
# 
#    classification = "C-CLOR2-MN-24-14"
# 
#    Problem data
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HIMMELBK'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['F'] = 142.22471
        v_['B1'] = 44.094
        v_['B2'] = 58.12
        v_['B3'] = 58.12
        v_['B4'] = 137.4
        v_['B5'] = 120.9
        v_['B6'] = 170.9
        v_['B7'] = 62.501
        v_['B8'] = 84.94
        v_['B9'] = 133.425
        v_['B10'] = 82.507
        v_['B11'] = 46.07
        v_['B12'] = 60.097
        v_['B13'] = 44.094
        v_['B14'] = 58.12
        v_['B15'] = 58.12
        v_['B16'] = 137.4
        v_['B17'] = 120.9
        v_['B18'] = 170.9
        v_['B19'] = 62.501
        v_['B20'] = 84.94
        v_['B21'] = 133.425
        v_['B22'] = 82.507
        v_['B23'] = 46.07
        v_['B24'] = 60.097
        v_['C1'] = 123.7
        v_['C2'] = 31.7
        v_['C3'] = 45.7
        v_['C4'] = 14.7
        v_['C5'] = 84.7
        v_['C6'] = 27.7
        v_['C7'] = 49.7
        v_['C8'] = 7.1
        v_['C9'] = 2.1
        v_['C10'] = 17.7
        v_['C11'] = 0.85
        v_['C12'] = 0.64
        v_['D1'] = 123.7
        v_['D2'] = 31.7
        v_['D3'] = 45.7
        v_['D4'] = 14.7
        v_['D5'] = 84.7
        v_['D6'] = 27.7
        v_['D7'] = 49.7
        v_['D8'] = 7.1
        v_['D9'] = 2.1
        v_['D10'] = 17.7
        v_['D11'] = 0.85
        v_['D12'] = 0.64
        v_['1'] = 1
        v_['12'] = 12
        v_['13'] = 13
        v_['24'] = 24
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for K in range(int(v_['1']),int(v_['24'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(K),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(K))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X1']])
        valA = np.append(valA,float(0.0693))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2']])
        valA = np.append(valA,float(0.0577))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X3']])
        valA = np.append(valA,float(0.05))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X4']])
        valA = np.append(valA,float(0.2))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X5']])
        valA = np.append(valA,float(0.26))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X6']])
        valA = np.append(valA,float(0.55))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X7']])
        valA = np.append(valA,float(0.06))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X8']])
        valA = np.append(valA,float(0.1))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X9']])
        valA = np.append(valA,float(0.12))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X10']])
        valA = np.append(valA,float(0.18))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X11']])
        valA = np.append(valA,float(0.1))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X12']])
        valA = np.append(valA,float(0.09))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X13']])
        valA = np.append(valA,float(0.0693))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X14']])
        valA = np.append(valA,float(0.0577))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X15']])
        valA = np.append(valA,float(0.05))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X16']])
        valA = np.append(valA,float(0.2))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X17']])
        valA = np.append(valA,float(0.26))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X18']])
        valA = np.append(valA,float(0.55))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X19']])
        valA = np.append(valA,float(0.06))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X20']])
        valA = np.append(valA,float(0.1))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X21']])
        valA = np.append(valA,float(0.12))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X22']])
        valA = np.append(valA,float(0.18))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X23']])
        valA = np.append(valA,float(0.1))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X24']])
        valA = np.append(valA,float(0.09))
        for I in range(int(v_['1']),int(v_['12'])+1):
            [ig,ig_,_] = s2mpj_ii('CA'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'CA'+str(I))
        for I in range(int(v_['1']),int(v_['24'])+1):
            [ig,ig_,_] = s2mpj_ii('CA13',ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'CA13')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(I)]])
            valA = np.append(valA,float(1.0))
        for I in range(int(v_['1']),int(v_['12'])+1):
            v_['I+12'] = 12+I
            v_['1/DI'] = 1.0/v_['D'+str(I)]
            [ig,ig_,_] = s2mpj_ii('CA14',ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'CA14')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(I)]])
            valA = np.append(valA,float(v_['1/DI']))
            v_['F/BI+12'] = v_['F']/v_['B'+str(int(v_['I+12']))]
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['I+12']))]])
            valA = np.append(valA,float(v_['F/BI+12']))
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
        self.gconst = arrset(self.gconst,ig_['CA13'],float(1.0))
        self.gconst = arrset(self.gconst,ig_['CA14'],float(1.671))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(0.04))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['1']),int(v_['12'])+1):
            v_['I+12'] = 12+I
            for J in range(int(v_['1']),int(v_['12'])+1):
                ename = 'E'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'en2PR')
                ielftype = arrset(ielftype,ie,iet_["en2PR"])
                vname = 'X'+str(int(v_['I+12']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.04))
                posev = np.where(elftv[ielftype[ie]]=='X')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'X'+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.04))
                posev = np.where(elftv[ielftype[ie]]=='Y')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
            for J in range(int(v_['13']),int(v_['24'])+1):
                ename = 'E'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'en2PR')
                ielftype = arrset(ielftype,ie,iet_["en2PR"])
                vname = 'X'+str(I)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.04))
                posev = np.where(elftv[ielftype[ie]]=='X')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'X'+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.04))
                posev = np.where(elftv[ielftype[ie]]=='Y')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['12'])+1):
            v_['I+12'] = 12+I
            for J in range(int(v_['1']),int(v_['12'])+1):
                v_['BI/BJ'] = v_['B'+str(I)]/v_['B'+str(J)]
                v_['40BI/BJ'] = 40.0*v_['BI/BJ']
                ig = ig_['CA'+str(I)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(v_['40BI/BJ']))
            for J in range(int(v_['13']),int(v_['24'])+1):
                v_['B+/BJ'] = v_['B'+str(int(v_['I+12']))]/v_['B'+str(J)]
                v_['CB+/BJ'] = v_['C'+str(I)]*v_['B+/BJ']
                v_['-CB+/BJ'] = -1.0*v_['CB+/BJ']
                ig = ig_['CA'+str(I)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(v_['-CB+/BJ']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN                0.0893344
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CLOR2-MN-24-14"
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

