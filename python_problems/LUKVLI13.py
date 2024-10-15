from s2mpjlib import *
class  LUKVLI13(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : LUKVLI13
#    *********
# 
#    Source: Problem 5.13, the chained modified HS48 problem, 
#    due to L. Luksan and J. Vlcek,
#    "Sparse and partially separable test problems for 
#    unconstrained and equality constrained optimization",
#    Technical Report 767, Inst. Computer Science, Academy of Sciences
#    of the Czech Republic, 182 07 Prague, Czech Republic, 1999
# 
#    Equality constraints changed to inequalities
# 
#    SIF input: Nick Gould, April 2001
# 
#    classification = "C-OQR2-AY-V-V"
# 
#    some useful parameters, including N, the number of variables.
# 
#           Alternative values for the SIF file parameters:
# IE N                   98             $-PARAMETER
# IE N                   998            $-PARAMETER
# IE N                   9998           $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'LUKVLI13'

    def __init__(self, *args): 
        import numpy as np
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(8);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
# IE N                   99998          $-PARAMETER
        v_['1'] = 1
        v_['2'] = 2
        v_['3'] = 3
        v_['4'] = 4
        v_['5'] = 5
        v_['N-2'] = -2+v_['N']
        v_['(N-2)/3'] = int(np.fix(v_['N-2']/v_['3']))
        v_['NC'] = v_['2']*v_['(N-2)/3']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.A       = lil_matrix((1000000,1000000))
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames      = np.array([])
        self.cnames = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['(N-2)/3'])+1):
            v_['I-1'] = -1+I
            v_['J'] = v_['3']*v_['I-1']
            v_['J+1'] = 1+v_['J']
            v_['J+2'] = 2+v_['J']
            v_['J+3'] = 3+v_['J']
            v_['J+4'] = 4+v_['J']
            v_['J+5'] = 5+v_['J']
            [ig,ig_,_] = s2mpj_ii('OBJ1'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(int(v_['J+1']))]
            self.A[ig,iv] = float(1.0)+self.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('OBJ2'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(int(v_['J+2']))]
            self.A[ig,iv] = float(1.0)+self.A[ig,iv]
            iv = ix_['X'+str(int(v_['J+3']))]
            self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('OBJ3'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(int(v_['J+4']))]
            self.A[ig,iv] = float(1.0)+self.A[ig,iv]
            iv = ix_['X'+str(int(v_['J+5']))]
            self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        for K in range(int(v_['1']),int(v_['NC'])+1,int(v_['2'])):
            v_['K+1'] = 1+K
            v_['K+2'] = 2+K
            v_['K+3'] = 3+K
            v_['K+4'] = 4+K
            [ig,ig_,_] = s2mpj_ii('C'+str(K),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'C'+str(K))
            iv = ix_['X'+str(K)]
            self.A[ig,iv] = float(1.0)+self.A[ig,iv]
            iv = ix_['X'+str(int(v_['K+2']))]
            self.A[ig,iv] = float(1.0)+self.A[ig,iv]
            iv = ix_['X'+str(int(v_['K+3']))]
            self.A[ig,iv] = float(1.0)+self.A[ig,iv]
            iv = ix_['X'+str(int(v_['K+4']))]
            self.A[ig,iv] = float(4.0)+self.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('C'+str(int(v_['K+1'])),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'C'+str(int(v_['K+1'])))
            iv = ix_['X'+str(int(v_['K+3']))]
            self.A[ig,iv] = float(-2.0)+self.A[ig,iv]
            iv = ix_['X'+str(int(v_['K+4']))]
            self.A[ig,iv] = float(-2.0)+self.A[ig,iv]
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
        for I in range(int(v_['1']),int(v_['(N-2)/3'])+1):
            self.gconst = arrset(self.gconst,ig_['OBJ1'+str(I)],float(1.0))
        for K in range(int(v_['1']),int(v_['NC'])+1,int(v_['2'])):
            v_['K+1'] = 1+K
            self.gconst = arrset(self.gconst,ig_['C'+str(K)],float(5.0))
            self.gconst = arrset(self.gconst,ig_['C'+str(int(v_['K+1']))],float(3.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        self.xlower = np.zeros((self.n,1))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        for I in range(int(v_['1']),int(v_['N'])+1,int(v_['3'])):
            self.x0[ix_['X'+str(I)]] = float(3.0)
        for I in range(int(v_['2']),int(v_['N'])+1,int(v_['3'])):
            self.x0[ix_['X'+str(I)]] = float(5.0)
        for I in range(int(v_['3']),int(v_['N'])+1,int(v_['3'])):
            self.x0[ix_['X'+str(I)]] = float(-3.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQR', iet_)
        elftv = loaset(elftv,it,0,'V')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for K in range(int(v_['1']),int(v_['NC'])+1,int(v_['2'])):
            v_['K+1'] = 1+K
            v_['K+2'] = 2+K
            ename = 'E'+str(K)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_["eSQR"])
            vname = 'X'+str(int(v_['K+1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'E'+str(int(v_['K+1']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_["eSQR"])
            ename = 'E'+str(int(v_['K+1']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(int(v_['K+2']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gL2',igt_)
        [it,igt_,_] = s2mpj_ii('gL4',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['(N-2)/3'])+1):
            ig = ig_['OBJ1'+str(I)]
            self.grftype = arrset(self.grftype,ig,'gL2')
            ig = ig_['OBJ2'+str(I)]
            self.grftype = arrset(self.grftype,ig,'gL2')
            ig = ig_['OBJ3'+str(I)]
            self.grftype = arrset(self.grftype,ig,'gL4')
        for K in range(int(v_['1']),int(v_['NC'])+1,int(v_['2'])):
            v_['K+1'] = 1+K
            ig = ig_['C'+str(K)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(K)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
            ig = ig_['C'+str(int(v_['K+1']))]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['K+1']))])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 0.0
#    Solution
# LO SOLTN               1.31224E+01
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.cupper[np.arange(self.nle)] = np.zeros((self.nle,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        self.A.resize(ngrp,self.n)
        self.A     = self.A.tocsr()
        sA1,sA2    = self.A.shape
        self.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass = "C-OQR2-AY-V-V"
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
        f_   = EV_[0]*EV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*EV_[0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0
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

    @staticmethod
    def gL4(self,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= GVAR_**4
        if nargout>1:
            g_ = 4.0*GVAR_**3
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 12.0*GVAR_**2
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

