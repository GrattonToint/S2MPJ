from s2mpjlib import *
class  PENALTY2(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : PENALTY2
#    --------
# 
#    The second penalty function
# 
#    This is a nonlinear least-squares problem with M=2*N groups.
#     Group 1 is linear.
#     Groups 2 to N use 2 nonlinear elements.
#     Groups N+1 to M-1 use 1 nonlinear element.
#     Group M uses N nonlinear elements.
#    The Hessian matrix is dense.
# 
#    Source:  Problem 24 in
#    J.J. More', B.S. Garbow and K.E. Hillstrom,
#    "Testing Unconstrained Optimization Software",
#    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
# 
#    See also Buckley#112 (p. 80)
# 
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "C-CSUR2-AN-V-0"
# 
#    Number of variables
# 
#           Alternative values for the SIF file parameters:
# IE N                   4              $-PARAMETER     original value
# IE N                   10             $-PARAMETER
# IE N                   50             $-PARAMETER
# IE N                   100            $-PARAMETER
# IE N                   200            $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'PENALTY2'

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
# IE N                   500            $-PARAMETER
# IE N                   1000           $-PARAMETER
        v_['A'] = 0.00001
        v_['B'] = 1.0
        v_['1'] = 1
        v_['2'] = 2
        v_['N+1'] = 1+v_['N']
        v_['M'] = v_['N']+v_['N']
        v_['M-1'] = -1+v_['M']
        v_['EM1/10'] = np.exp(-0.1)
        v_['1/A'] = 1.0/v_['A']
        v_['1/B'] = 1.0/v_['B']
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
        [ig,ig_,_] = s2mpj_ii('G'+str(int(v_['1'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X'+str(int(v_['1']))]])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('G'+str(int(v_['1'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        self.gscale = arrset(self.gscale,ig,float(v_['1/B']))
        for I in range(int(v_['2']),int(v_['M-1'])+1):
            [ig,ig_,_] = s2mpj_ii('G'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            self.gscale = arrset(self.gscale,ig,float(v_['1/A']))
        [ig,ig_,_] = s2mpj_ii('G'+str(int(v_['M'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        self.gscale = arrset(self.gscale,ig,float(v_['1/B']))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        selfnob      = ngrp
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        self.gconst = arrset(self.gconst,ig_['G1'],float(0.2))
        for I in range(int(v_['2']),int(v_['N'])+1):
            v_['I-1'] = -1+I
            v_['RI'] = float(I)
            v_['RI-1'] = float(v_['I-1'])
            v_['I/10'] = 0.1*v_['RI']
            v_['I-1/10'] = 0.1*v_['RI-1']
            v_['EI/10'] = np.exp(v_['I/10'])
            v_['EI-1/10'] = np.exp(v_['I-1/10'])
            v_['YI'] = v_['EI/10']+v_['EI-1/10']
            self.gconst = arrset(self.gconst,ig_['G'+str(I)],float(v_['YI']))
        for I in range(int(v_['N+1']),int(v_['M-1'])+1):
            self.gconst = arrset(self.gconst,ig_['G'+str(I)],float(v_['EM1/10']))
        self.gconst = arrset(self.gconst,ig_['G'+str(int(v_['M']))],float(1.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(0.5))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eE10', iet_)
        elftv = loaset(elftv,it,0,'V')
        [it,iet_,_] = s2mpj_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'V')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['2']),int(v_['N'])+1):
            v_['I-1'] = -1+I
            ename = 'A'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eE10')
            ielftype = arrset(ielftype,ie,iet_["eE10"])
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.5))
            posev = np.where(elftv[ielftype[ie]]=='V')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eE10')
            ielftype = arrset(ielftype,ie,iet_["eE10"])
            vname = 'X'+str(int(v_['I-1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.5))
            posev = np.where(elftv[ielftype[ie]]=='V')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        for I in range(int(v_['N+1']),int(v_['M-1'])+1):
            v_['-N'] = -1*v_['N']
            v_['I-N'] = I+v_['-N']
            v_['I-N+1'] = 1+v_['I-N']
            ename = 'C'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eE10')
            ielftype = arrset(ielftype,ie,iet_["eE10"])
            vname = 'X'+str(int(v_['I-N+1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.5))
            posev = np.where(elftv[ielftype[ie]]=='V')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        for J in range(int(v_['1']),int(v_['N'])+1):
            ename = 'D'+str(J)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSQ')
            ielftype = arrset(ielftype,ie,iet_["eSQ"])
            vname = 'X'+str(J)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.5))
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
        ig = ig_['G'+str(int(v_['1']))]
        self.grftype = arrset(self.grftype,ig,'gL2')
        for I in range(int(v_['2']),int(v_['N'])+1):
            ig = ig_['G'+str(I)]
            self.grftype = arrset(self.grftype,ig,'gL2')
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['A'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = posel+1
            self.grelt = loaset(self.grelt,ig,posel,ie_['B'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel, 1.)
        for I in range(int(v_['N+1']),int(v_['M-1'])+1):
            ig = ig_['G'+str(I)]
            self.grftype = arrset(self.grftype,ig,'gL2')
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['C'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['G'+str(int(v_['M']))]
        self.grftype = arrset(self.grftype,ig,'gL2')
        for J in range(int(v_['1']),int(v_['N'])+1):
            v_['-J'] = -1*J
            v_['N-J'] = v_['N']+v_['-J']
            v_['N-J+1'] = 1+v_['N-J']
            v_['WI'] = float(v_['N-J+1'])
            ig = ig_['G'+str(int(v_['M']))]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['D'+str(J)])
            self.grelw = loaset(self.grelw,ig,posel,float(v_['WI']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        self.objlower = 0.0
#    Solution
# LO SOLTN(4)            9.37629D-6
# LO SOLTN(10)           2.93660D-4
# LO SOLTN(50)           4.29609813
# LO SOLTN(100)          97096.0840
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass   = "C-CSUR2-AN-V-0"
        self.objderlvl = 2

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eE10(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        EXPA = np.exp(0.1*EV_[0,0])
        f_   = EXPA
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 0.1*EXPA
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 0.01*EXPA
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eSQ(self, nargout,*args):

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
            g_[0] = EV_[0,0]+EV_[0,0]
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

