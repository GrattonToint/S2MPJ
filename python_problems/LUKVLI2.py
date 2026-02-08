from s2mpjlib import *
class  LUKVLI2(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : LUKVLI2
#    *********
# 
#    Source: Problem 5.2, the chained Wood function with 
#    Broyden banded constraints, due to L. Luksan and J. Vlcek,
#    "Sparse and partially separable test problems for 
#    unconstrained and equality constrained optimization",
#    Technical Report 767, Inst. Computer Science, Academy of Sciences
#    of the Czech Republic, 182 07 Prague, Czech Republic, 1999
# 
#    Equality constraints changed to inequalities
# 
#    SIF input: Nick Gould, April 2001
# 
#    classification = "C-COOR2-AY-V-V"
# 
#    some useful parameters, including N, the number of variables.
# 
#           Alternative values for the SIF file parameters:
# IE N                   100            $-PARAMETER
# IE N                   1000           $-PARAMETER
# IE N                   10000          $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'LUKVLI2'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(10);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
# IE N                   100000         $-PARAMETER
        v_['0'] = 0
        v_['1'] = 1
        v_['2'] = 2
        v_['6'] = 6
        v_['7'] = 7
        v_['N/2'] = int(np.fix(v_['N']/v_['2']))
        v_['N/2-1'] = -1+v_['N/2']
        v_['N-1'] = -1+v_['N']
        v_['N-2'] = -2+v_['N']
        v_['1.0'] = 1.0
        v_['90.0'] = 90.0
        v_['1/90'] = v_['1.0']/v_['90.0']
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
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['1']),int(v_['N/2-1'])+1):
            v_['2I'] = 2*I
            v_['2I-1'] = -1+v_['2I']
            v_['2I+1'] = 1+v_['2I']
            v_['2I+2'] = 2+v_['2I']
            [ig,ig_,_] = s2mpj_ii('A'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['2I']))]])
            valA = np.append(valA,float(-1.0))
            self.gscale = arrset(self.gscale,ig,float(0.01))
            [ig,ig_,_] = s2mpj_ii('B'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['2I-1']))]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('C'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['2I+2']))]])
            valA = np.append(valA,float(-1.0))
            self.gscale = arrset(self.gscale,ig,float(v_['1/90']))
            [ig,ig_,_] = s2mpj_ii('D'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['2I+1']))]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('E'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['2I']))]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['2I+2']))]])
            valA = np.append(valA,float(1.0))
            self.gscale = arrset(self.gscale,ig,float(0.1))
            [ig,ig_,_] = s2mpj_ii('F'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['2I']))]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['2I-1']))]])
            valA = np.append(valA,float(-1.0))
            self.gscale = arrset(self.gscale,ig,float(10.0))
        for K in range(int(v_['6']),int(v_['N-2'])+1):
            v_['K+1'] = 1+K
            v_['K-5'] = -5+K
            [ig,ig_,_] = s2mpj_ii('C'+str(K),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'C'+str(K))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(K)]])
            valA = np.append(valA,float(2.0))
            for I in range(int(v_['K-5']),int(v_['K+1'])+1):
                [ig,ig_,_] = s2mpj_ii('C'+str(K),ig_)
                gtype = arrset(gtype,ig,'<=')
                cnames = arrset(cnames,ig,'C'+str(K))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(I)]])
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
        for I in range(int(v_['1']),int(v_['N/2-1'])+1):
            self.gconst = arrset(self.gconst,ig_['B'+str(I)],float(1.0))
            self.gconst = arrset(self.gconst,ig_['D'+str(I)],float(1.0))
            self.gconst = arrset(self.gconst,ig_['E'+str(I)],float(2.0))
        for K in range(int(v_['6']),int(v_['N-2'])+1):
            self.gconst = arrset(self.gconst,ig_['C'+str(K)],float(-1.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        for I in range(int(v_['1']),int(v_['N'])+1,int(v_['2'])):
            self.x0[ix_['X'+str(I)]] = float(-2.0)
        for I in range(int(v_['2']),int(v_['N'])+1,int(v_['2'])):
            self.x0[ix_['X'+str(I)]] = float(1.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQR', iet_)
        elftv = loaset(elftv,it,0,'V')
        [it,iet_,_] = s2mpj_ii( 'eCUBE', iet_)
        elftv = loaset(elftv,it,0,'V')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['1']),int(v_['N/2-1'])+1):
            v_['2I'] = 2*I
            v_['2I-1'] = -1+v_['2I']
            v_['2I+1'] = 1+v_['2I']
            ename = 'A'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_["eSQR"])
            vname = 'X'+str(int(v_['2I-1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'C'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_["eSQR"])
            vname = 'X'+str(int(v_['2I+1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        for K in range(int(v_['6']),int(v_['N-2'])+1):
            ename = 'U'+str(K)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eCUBE')
            ielftype = arrset(ielftype,ie,iet_["eCUBE"])
            vname = 'X'+str(K)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        for K in range(int(v_['1']),int(v_['N-1'])+1):
            ename = 'S'+str(K)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_["eSQR"])
            vname = 'X'+str(K)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gL2',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['N/2-1'])+1):
            ig = ig_['A'+str(I)]
            self.grftype = arrset(self.grftype,ig,'gL2')
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['A'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
            ig = ig_['B'+str(I)]
            self.grftype = arrset(self.grftype,ig,'gL2')
            ig = ig_['C'+str(I)]
            self.grftype = arrset(self.grftype,ig,'gL2')
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['C'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
            ig = ig_['D'+str(I)]
            self.grftype = arrset(self.grftype,ig,'gL2')
            ig = ig_['E'+str(I)]
            self.grftype = arrset(self.grftype,ig,'gL2')
            ig = ig_['F'+str(I)]
            self.grftype = arrset(self.grftype,ig,'gL2')
        for K in range(int(v_['6']),int(v_['N-2'])+1):
            v_['K+1'] = 1+K
            v_['K-5'] = -5+K
            ig = ig_['C'+str(K)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['U'+str(K)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(5.0))
            for I in range(int(v_['K-5']),int(v_['K+1'])+1):
                ig = ig_['C'+str(K)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['S'+str(I)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 0.0
#    Solution
# LO SOLTN               2.39390E+04
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
        self.pbclass   = "C-COOR2-AY-V-V"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

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
            g_[0] = 2.0*EV_[0,0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eCUBE(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0,0]**3
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 3.0*EV_[0,0]**2
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 6.0*EV_[0,0]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gL2(self,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= GVAR_*GVAR_
        if nargout>1:
            g_ = GVAR_+GVAR_
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

