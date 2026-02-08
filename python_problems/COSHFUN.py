from s2mpjlib import *
class  COSHFUN(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : COSHFUN
#    *********
# 
#    A nonlinear minmax problem.
# 
#    Source:
#    K. Jonasson and K. Madsen,
#    "Corrected sequential linear programming for sparse
#    minimax optimization", Technical report, Institute for Numerical
#    Analysis, Technical U. of Denmark.
# 
#    SIF input: Nick Gould, October 1992.
# 
#    classification = "C-CLOR2-AN-V-V"
# 
#   the number of functions
# 
#           Alternative values for the SIF file parameters:
# IE M                   3              $-PARAMETER
# IE M                   8              $-PARAMETER
# IE M                   14             $-PARAMETER
# IE M                   20             $-PARAMETER     original value
# IE M                   200            $-PARAMETER
# IE M                   2000           $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'COSHFUN'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['M'] = int(8);  #  SIF file default value
        else:
            v_['M'] = int(args[0])
        v_['N'] = 3*v_['M']
        v_['N-3'] = -3+v_['N']
        v_['N-5'] = -5+v_['N']
        v_['0'] = 0
        v_['1'] = 1
        v_['3'] = 3
        v_['6'] = 6
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
        [iv,ix_,_] = s2mpj_ii('F',ix_)
        self.xnames=arrset(self.xnames,iv,'F')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['F']])
        valA = np.append(valA,float(1.0))
        for I in range(int(v_['1']),int(v_['M'])+1):
            [ig,ig_,_] = s2mpj_ii('C'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'C'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['F']])
            valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('C1',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'C1')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X3']])
        valA = np.append(valA,float(-2.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X6']])
        valA = np.append(valA,float(-1.0))
        for I in range(int(v_['6']),int(v_['N-3'])+1,int(v_['3'])):
            v_['I-5'] = -5+I
            v_['I+3'] = 3+I
            v_['I/3'] = int(np.fix(I/v_['3']))
            [ig,ig_,_] = s2mpj_ii('C'+str(int(v_['I/3'])),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'C'+str(int(v_['I/3'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['I-5']))]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(I)]])
            valA = np.append(valA,float(-2.0))
            [ig,ig_,_] = s2mpj_ii('C'+str(int(v_['I/3'])),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'C'+str(int(v_['I/3'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['I+3']))]])
            valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('C'+str(int(v_['M'])),ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'C'+str(int(v_['M'])))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X'+str(int(v_['N-5']))]])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X'+str(int(v_['N']))]])
        valA = np.append(valA,float(-2.0))
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
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQR', iet_)
        elftv = loaset(elftv,it,0,'X')
        [it,iet_,_] = s2mpj_ii( 'eCOSH', iet_)
        elftv = loaset(elftv,it,0,'X')
        [it,iet_,_] = s2mpj_ii( 'ePROD', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['3']),int(v_['N'])+1,int(v_['3'])):
            v_['I-1'] = -1+I
            v_['I-2'] = -2+I
            v_['I/3'] = int(np.fix(I/v_['3']))
            ename = 'SQR'+str(int(v_['I/3']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_["eSQR"])
            ename = 'SQR'+str(int(v_['I/3']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.x0 = np.zeros((self.n,1))
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'COSH'+str(int(v_['I/3']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eCOSH')
            ielftype = arrset(ielftype,ie,iet_["eCOSH"])
            ename = 'COSH'+str(int(v_['I/3']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(int(v_['I-1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'PROD'+str(int(v_['I/3']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'ePROD')
            ielftype = arrset(ielftype,ie,iet_["ePROD"])
            ename = 'PROD'+str(int(v_['I/3']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(int(v_['I-2']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'PROD'+str(int(v_['I/3']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='Y')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['M'])+1):
            ig = ig_['C'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['SQR'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = posel+1
            self.grelt = loaset(self.grelt,ig,posel,ie_['COSH'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel, 1.)
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['PROD'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.cupper[np.arange(self.nle)] = np.zeros((self.nle,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CLOR2-AN-V-V"
        self.x0        = np.zeros((self.n,1))
        self.objderlvl = 2
        self.conderlvl = [2]


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSQR(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0,0]*EV_[0,0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0e+0*EV_[0,0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0e+0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eCOSH(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        COSHX = np.cosh(EV_[0,0])
        f_   = COSHX
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = np.sinh(EV_[0,0])
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = COSHX
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ePROD(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = 2.0e+0*EV_[0,0]*EV_[0,0]*EV_[1,0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 4.0e+0*EV_[0,0]*EV_[1,0]
            g_[1] = 2.0e+0*EV_[0,0]*EV_[0,0]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 4.0e+0*EV_[1,0]
                H_[0,1] = 4.0e+0*EV_[0,0]
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

