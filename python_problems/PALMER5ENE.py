from s2mpjlib import *
class  PALMER5ENE(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : PALMER5ENE
#    *********
# 
#    A nonlinear least squares problem
#    arising from chemical kinetics.
# 
#    model: H-N=C=Se TZVP + MP2
#    fitting Y to A0 T_0 + A2 T_2 + A4 T_4 + A6 T_6 + A8 T_8 +
#                 A10 T_10 + A12 T_12 + A14 T_14
#                 + L * EXP( -K X**2 )
#    where T_i is the i-th (shifted) Chebyshev polynomial
# 
#    Source:
#    M.  Palmer, Edinburgh, private communication.
# 
#    SIF input: Nick Gould, 1992.
#    Bound-constrained nonlinear equations version: Nick Gould, June 2019.
# 
#    classification = "C-CNOR2-RN-8-12"
# 
#    Number of data points
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 25 XI 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'PALMER5ENE'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['M'] = 23
        v_['1'] = 1
        v_['2'] = 2
        v_['12'] = 12
        v_['14'] = 14
        v_['X12'] = 0.000000
        v_['X13'] = 1.570796
        v_['X14'] = 1.396263
        v_['X15'] = 1.308997
        v_['X16'] = 1.221730
        v_['X17'] = 1.125835
        v_['X18'] = 1.047198
        v_['X19'] = 0.872665
        v_['X20'] = 0.698132
        v_['X21'] = 0.523599
        v_['X22'] = 0.349066
        v_['X23'] = 0.174533
        v_['B'] = v_['X13']
        v_['A'] = -1.0e+0*v_['B']
        v_['DIFF'] = 2.0e+0*v_['B']
        v_['Y12'] = 83.57418
        v_['Y13'] = 81.007654
        v_['Y14'] = 18.983286
        v_['Y15'] = 8.051067
        v_['Y16'] = 2.044762
        v_['Y17'] = 0.000000
        v_['Y18'] = 1.170451
        v_['Y19'] = 10.479881
        v_['Y20'] = 25.785001
        v_['Y21'] = 44.126844
        v_['Y22'] = 62.822177
        v_['Y23'] = 77.719674
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [iv,ix_,_] = s2mpj_ii('A0',ix_)
        self.xnames=arrset(self.xnames,iv,'A0')
        [iv,ix_,_] = s2mpj_ii('A2',ix_)
        self.xnames=arrset(self.xnames,iv,'A2')
        [iv,ix_,_] = s2mpj_ii('A4',ix_)
        self.xnames=arrset(self.xnames,iv,'A4')
        [iv,ix_,_] = s2mpj_ii('A6',ix_)
        self.xnames=arrset(self.xnames,iv,'A6')
        [iv,ix_,_] = s2mpj_ii('A8',ix_)
        self.xnames=arrset(self.xnames,iv,'A8')
        [iv,ix_,_] = s2mpj_ii('A10',ix_)
        self.xnames=arrset(self.xnames,iv,'A10')
        [iv,ix_,_] = s2mpj_ii('K',ix_)
        self.xnames=arrset(self.xnames,iv,'K')
        [iv,ix_,_] = s2mpj_ii('L',ix_)
        self.xnames=arrset(self.xnames,iv,'L')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['12']),int(v_['M'])+1):
            v_['T0'] = 1.0e+0
            v_['Y'] = 2.0e+0*v_['X'+str(I)]
            v_['Y'] = v_['Y']-v_['A']
            v_['Y'] = v_['Y']-v_['B']
            v_['Y'] = v_['Y']/v_['DIFF']
            v_['T1'] = v_['Y']
            v_['2Y'] = 2.0e+0*v_['Y']
            for J in range(int(v_['2']),int(v_['14'])+1):
                v_['J-1'] = -1+J
                v_['J-2'] = -2+J
                v_['T'+str(J)] = v_['2Y']*v_['T'+str(int(v_['J-1']))]
                v_['T'+str(J)] = v_['T'+str(J)]-v_['T'+str(int(v_['J-2']))]
            [ig,ig_,_] = s2mpj_ii('O'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'O'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['A0']])
            valA = np.append(valA,float(v_['T0']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['A2']])
            valA = np.append(valA,float(v_['T2']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['A4']])
            valA = np.append(valA,float(v_['T4']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['A6']])
            valA = np.append(valA,float(v_['T6']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['A8']])
            valA = np.append(valA,float(v_['T8']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['A10']])
            valA = np.append(valA,float(v_['T10']))
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
        for I in range(int(v_['12']),int(v_['M'])+1):
            self.gconst = arrset(self.gconst,ig_['O'+str(I)],float(v_['Y'+str(I)]))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        self.xlower[ix_['A0']] = -float('Inf')
        self.xupper[ix_['A0']] = +float('Inf')
        self.xlower[ix_['A2']] = -float('Inf')
        self.xupper[ix_['A2']] = +float('Inf')
        self.xlower[ix_['A4']] = -float('Inf')
        self.xupper[ix_['A4']] = +float('Inf')
        self.xlower[ix_['A6']] = -float('Inf')
        self.xupper[ix_['A6']] = +float('Inf')
        self.xlower[ix_['A8']] = -float('Inf')
        self.xupper[ix_['A8']] = +float('Inf')
        self.xlower[ix_['A10']] = -float('Inf')
        self.xupper[ix_['A10']] = +float('Inf')
        self.xlower[ix_['L']] = -float('Inf')
        self.xupper[ix_['L']] = +float('Inf')
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(1.0))
        self.x0[ix_['A0']] = float(1.9264e+01)
        self.x0[ix_['A2']] = float(-1.7302)
        self.x0[ix_['A4']] = float(4.0794e+01)
        self.x0[ix_['A6']] = float(8.3021e-01)
        self.x0[ix_['A8']] = float(3.7090)
        self.x0[ix_['A10']] = float(-1.7723e-01)
        self.x0[ix_['K']] = float(10.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'ePROD', iet_)
        elftv = loaset(elftv,it,0,'K')
        elftv = loaset(elftv,it,1,'L')
        elftp = []
        elftp = loaset(elftp,it,0,'XSQR')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['12']),int(v_['M'])+1):
            v_['XSQR'] = v_['X'+str(I)]*v_['X'+str(I)]
            ename = 'E'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'ePROD')
            ielftype = arrset(ielftype,ie,iet_["ePROD"])
            vname = 'K'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(1.0))
            posev = np.where(elftv[ielftype[ie]]=='K')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'L'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(1.0))
            posev = np.where(elftv[ielftype[ie]]=='L')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='XSQR')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['XSQR']))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['12']),int(v_['M'])+1):
            ig = ig_['O'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
# LO PALMER5E                0.0
#    Solution
# LO SOLTN              1.48003482D-04
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
        self.pbclass   = "C-CNOR2-RN-8-12"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def ePROD(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        EXPON = np.exp(-EV_[0]*self.elpar[iel_][0])
        f_   = EV_[1]*EXPON
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -self.elpar[iel_][0]*EV_[1]*EXPON
            g_[1] = EXPON
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = self.elpar[iel_][0]*self.elpar[iel_][0]*EV_[1]*EXPON
                H_[0,1] = -self.elpar[iel_][0]*EXPON
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

