from s2mpjlib import *
class  PALMER1BNE(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : PALMER1BNE
#    *********
# 
#    A nonlinear least squares problem with bounds
#    arising from chemical kinetics.
# 
#    model: H-N=N=N TZVP+MP2
#    fitting Y to A2 X**2 + A4 X**4
#                 + B / ( C + X**2 ), B, C nonnegative.
# 
#    Source:
#    M. Palmer, Edinburgh, private communication.
# 
#    SIF input: Nick Gould, 1990.
#    Bound-constrained nonlinear equations version: Nick Gould, June 2019.
# 
#    classification = "C-CNOR2-RN-4-35"
# 
#    Number of data points
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'PALMER1BNE'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['M'] = 35
        v_['1'] = 1
        v_['X1'] = -1.788963
        v_['X2'] = -1.745329
        v_['X3'] = -1.658063
        v_['X4'] = -1.570796
        v_['X5'] = -1.483530
        v_['X6'] = -1.396263
        v_['X7'] = -1.308997
        v_['X8'] = -1.218612
        v_['X9'] = -1.134464
        v_['X10'] = -1.047198
        v_['X11'] = -0.872665
        v_['X12'] = -0.698132
        v_['X13'] = -0.523599
        v_['X14'] = -0.349066
        v_['X15'] = -0.174533
        v_['X16'] = 0.0000000
        v_['X17'] = 1.788963
        v_['X18'] = 1.745329
        v_['X19'] = 1.658063
        v_['X20'] = 1.570796
        v_['X21'] = 1.483530
        v_['X22'] = 1.396263
        v_['X23'] = 1.308997
        v_['X24'] = 1.218612
        v_['X25'] = 1.134464
        v_['X26'] = 1.047198
        v_['X27'] = 0.872665
        v_['X28'] = 0.698132
        v_['X29'] = 0.523599
        v_['X30'] = 0.349066
        v_['X31'] = 0.174533
        v_['X32'] = -1.8762289
        v_['X33'] = -1.8325957
        v_['X34'] = 1.8762289
        v_['X35'] = 1.8325957
        v_['Y1'] = 78.596218
        v_['Y2'] = 65.77963
        v_['Y3'] = 43.96947
        v_['Y4'] = 27.038816
        v_['Y5'] = 14.6126
        v_['Y6'] = 6.2614
        v_['Y7'] = 1.538330
        v_['Y8'] = 0.000000
        v_['Y9'] = 1.188045
        v_['Y10'] = 4.6841
        v_['Y11'] = 16.9321
        v_['Y12'] = 33.6988
        v_['Y13'] = 52.3664
        v_['Y14'] = 70.1630
        v_['Y15'] = 83.4221
        v_['Y16'] = 88.3995
        v_['Y17'] = 78.596218
        v_['Y18'] = 65.77963
        v_['Y19'] = 43.96947
        v_['Y20'] = 27.038816
        v_['Y21'] = 14.6126
        v_['Y22'] = 6.2614
        v_['Y23'] = 1.538330
        v_['Y24'] = 0.000000
        v_['Y25'] = 1.188045
        v_['Y26'] = 4.6841
        v_['Y27'] = 16.9321
        v_['Y28'] = 33.6988
        v_['Y29'] = 52.3664
        v_['Y30'] = 70.1630
        v_['Y31'] = 83.4221
        v_['Y32'] = 108.18086
        v_['Y33'] = 92.733676
        v_['Y34'] = 108.18086
        v_['Y35'] = 92.733676
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [iv,ix_,_] = s2mpj_ii('A2',ix_)
        self.xnames=arrset(self.xnames,iv,'A2')
        [iv,ix_,_] = s2mpj_ii('A4',ix_)
        self.xnames=arrset(self.xnames,iv,'A4')
        [iv,ix_,_] = s2mpj_ii('B',ix_)
        self.xnames=arrset(self.xnames,iv,'B')
        [iv,ix_,_] = s2mpj_ii('C',ix_)
        self.xnames=arrset(self.xnames,iv,'C')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['1']),int(v_['M'])+1):
            v_['XSQR'] = v_['X'+str(I)]*v_['X'+str(I)]
            v_['XQUART'] = v_['XSQR']*v_['XSQR']
            [ig,ig_,_] = s2mpj_ii('O'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'O'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['A2']])
            valA = np.append(valA,float(v_['XSQR']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['A4']])
            valA = np.append(valA,float(v_['XQUART']))
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
        for I in range(int(v_['1']),int(v_['M'])+1):
            self.gconst = arrset(self.gconst,ig_['O'+str(I)],float(v_['Y'+str(I)]))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        self.xlower[ix_['A2']] = -float('Inf')
        self.xupper[ix_['A2']] = +float('Inf')
        self.xlower[ix_['A4']] = -float('Inf')
        self.xupper[ix_['A4']] = +float('Inf')
        self.xlower[ix_['B']] = 0.00001
        self.xlower[ix_['C']] = 0.00001
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(1.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eQUOT', iet_)
        elftv = loaset(elftv,it,0,'B')
        elftv = loaset(elftv,it,1,'C')
        elftp = []
        elftp = loaset(elftp,it,0,'XSQR')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['1']),int(v_['M'])+1):
            v_['XSQR'] = v_['X'+str(I)]*v_['X'+str(I)]
            ename = 'E'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eQUOT')
            ielftype = arrset(ielftype,ie,iet_["eQUOT"])
            vname = 'B'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(1.0))
            posev = np.where(elftv[ielftype[ie]]=='B')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'C'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(1.0))
            posev = np.where(elftv[ielftype[ie]]=='C')[0]
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
        for I in range(int(v_['1']),int(v_['M'])+1):
            ig = ig_['O'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
# LO PALMER1B               0.0
#    Solution
# LO SOLTN               3.44734948
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
        self.pbclass   = "C-CNOR2-RN-4-35"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eQUOT(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        DENOM = 1.0/(EV_[1,0]+self.elpar[iel_][0])
        f_   = EV_[0,0]*DENOM
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = DENOM
            g_[1] = -EV_[0,0]*DENOM*DENOM
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = -DENOM*DENOM
                H_[1,0] = H_[0,1]
                H_[1,1] = 2.0*EV_[0,0]*DENOM**3
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

