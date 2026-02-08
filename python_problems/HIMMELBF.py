from s2mpjlib import *
class  HIMMELBF(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HIMMELBF
#    *********
# 
#    A 4 variables data fitting problems by Himmelblau.
# 
#    Source: problem 32 in
#    D.H. Himmelblau,
#    "Applied nonlinear programming",
#    McGraw-Hill, New-York, 1972.
# 
#    See Buckley#76 (p. 66)
# 
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "C-CSUR2-AN-4-0"
# 
#    Problem data
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HIMMELBF'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['A1'] = 0.0
        v_['A2'] = 0.000428
        v_['A3'] = 0.001000
        v_['A4'] = 0.001610
        v_['A5'] = 0.002090
        v_['A6'] = 0.003480
        v_['A7'] = 0.005250
        v_['B1'] = 7.391
        v_['B2'] = 11.18
        v_['B3'] = 16.44
        v_['B4'] = 16.20
        v_['B5'] = 22.20
        v_['B6'] = 24.02
        v_['B7'] = 31.32
        v_['1'] = 1
        v_['7'] = 7
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [iv,ix_,_] = s2mpj_ii('X1',ix_)
        self.xnames=arrset(self.xnames,iv,'X1')
        [iv,ix_,_] = s2mpj_ii('X2',ix_)
        self.xnames=arrset(self.xnames,iv,'X2')
        [iv,ix_,_] = s2mpj_ii('X3',ix_)
        self.xnames=arrset(self.xnames,iv,'X3')
        [iv,ix_,_] = s2mpj_ii('X4',ix_)
        self.xnames=arrset(self.xnames,iv,'X4')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['1']),int(v_['7'])+1):
            [ig,ig_,_] = s2mpj_ii('G'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            self.gscale = arrset(self.gscale,ig,float(0.0001))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        selfnob      = ngrp
        #%%%%%%%%%%%%%%%%%%  CONSTANTS %%%%%%%%%%%%%%%%%%%
        self.gconst = np.full((ngrp,1),1.0)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.x0[ix_['X1']] = float(2.7)
        self.x0[ix_['X2']] = float(90.0)
        self.x0[ix_['X3']] = float(1500.0)
        self.x0[ix_['X4']] = float(10.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eHF', iet_)
        elftv = loaset(elftv,it,0,'XA')
        elftv = loaset(elftv,it,1,'XB')
        elftv = loaset(elftv,it,2,'XC')
        elftv = loaset(elftv,it,3,'XD')
        elftp = []
        elftp = loaset(elftp,it,0,'A')
        elftp = loaset(elftp,it,1,'B')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['1']),int(v_['7'])+1):
            ename = 'E'+str(I)
            [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
            if newelt:
                self.elftype = arrset(self.elftype,ie,'eHF')
                ielftype = arrset(ielftype,ie,iet_['eHF'])
            vname = 'X1'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='XA')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='XB')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='XC')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X4'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='XD')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='A')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['A'+str(I)]))
            posep = np.where(elftp[ielftype[ie]]=='B')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['B'+str(I)]))
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
        for ig in range(0,ngrp):
            self.grftype = arrset(self.grftype,ig,'gL2')
        for I in range(int(v_['1']),int(v_['7'])+1):
            ig = ig_['G'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               318.572
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass   = "C-CSUR2-AN-4-0"
        self.objderlvl = 2

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eHF(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U  = (
              EV_[0,0]*EV_[0,0]+self.elpar[iel_][0]*EV_[1,0]*EV_[1,0]+self.elpar[iel_][0]*self.elpar[iel_][0]*EV_[2,0]*EV_[2,0])
        V = self.elpar[iel_][1]*(1.0+self.elpar[iel_][0]*EV_[3,0]*EV_[3,0])
        V2 = V*V
        AB = self.elpar[iel_][0]*self.elpar[iel_][1]
        A2 = self.elpar[iel_][0]*self.elpar[iel_][0]
        T = -4.0*AB/V2
        f_   = U/V
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*EV_[0,0]/V
            g_[1] = 2.0*self.elpar[iel_][0]*EV_[1,0]/V
            g_[2] = 2.0*A2*EV_[2,0]/V
            g_[3] = -2.0*AB*EV_[3,0]*U/V2
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,0] = 2.0/V
                H_[0,3] = T*EV_[3,0]*EV_[0,0]
                H_[3,0] = H_[0,3]
                H_[1,1] = 2.0*self.elpar[iel_][0]/V
                H_[1,3] = T*self.elpar[iel_][0]*EV_[3,0]*EV_[1,0]
                H_[3,1] = H_[1,3]
                H_[2,2] = 2.0*A2/V
                H_[2,3] = T*self.elpar[iel_][0]*EV_[3,0]*EV_[2,0]
                H_[3,2] = H_[2,3]
                H_[3,3] = -2.0*AB*U/V2+8.0*(AB*EV_[3,0])**2*U/(V2*V)
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

