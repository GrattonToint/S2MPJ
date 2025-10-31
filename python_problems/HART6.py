from s2mpjlib import *
class  HART6(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HART6
#    *********
# 
#    Source: Hartman problem 6 in
#    L. C. W. Dixon and G. P. Szego (Eds.)
#    Towards Global Optimization
#    North Holland, 1975.
#    Paper 9, page 163.
# 
#    SIF input: A.R. Conn May 1995
# 
#    classification = "C-COBR2-AN-6-0"
# 
#    Number of variables - constraints
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HART6'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 6
        v_['1'] = 1
        v_['ONE'] = 1
        v_['NN'] = 6
        v_['L'] = 4
        v_['C1'] = 1.0
        v_['C2'] = 1.2
        v_['C3'] = 3.0
        v_['C4'] = 3.2
        v_['A1,1'] = 10.0
        v_['A2,1'] = 0.05
        v_['A3,1'] = 3.0
        v_['A4,1'] = 17.0
        v_['A1,2'] = 0.05
        v_['A2,2'] = 10.0
        v_['A3,2'] = 3.5
        v_['A4,2'] = 8.0
        v_['A1,3'] = 17.0
        v_['A2,3'] = 17.0
        v_['A3,3'] = 1.7
        v_['A4,3'] = 0.05
        v_['A1,4'] = 3.5
        v_['A2,4'] = 0.1
        v_['A3,4'] = 10.0
        v_['A4,4'] = 10.0
        v_['A1,5'] = 1.7
        v_['A2,5'] = 8.0
        v_['A3,5'] = 17.0
        v_['A4,5'] = 0.1
        v_['A1,6'] = 8.0
        v_['A2,6'] = 14.0
        v_['A3,6'] = 8.0
        v_['A4,6'] = 14.0
        v_['P1,1'] = 0.1312
        v_['P2,1'] = 0.2329
        v_['P3,1'] = 0.2348
        v_['P4,1'] = 0.4047
        v_['P1,2'] = 0.1696
        v_['P2,2'] = 0.4135
        v_['P3,2'] = 0.1451
        v_['P4,2'] = 0.8828
        v_['P1,3'] = 0.5569
        v_['P2,3'] = 0.8307
        v_['P3,3'] = 0.3522
        v_['P4,3'] = 0.8732
        v_['P1,4'] = 0.0124
        v_['P2,4'] = 0.3736
        v_['P3,4'] = 0.2883
        v_['P4,4'] = 0.5743
        v_['P1,5'] = 0.8283
        v_['P2,5'] = 0.1004
        v_['P3,5'] = 0.3047
        v_['P4,5'] = 0.1091
        v_['P1,6'] = 0.5886
        v_['P2,6'] = 0.9991
        v_['P3,6'] = 0.6650
        v_['P4,6'] = 0.0381
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
        for I in range(int(v_['1']),int(v_['L'])+1):
            [ig,ig_,_] = s2mpj_ii('OBJ'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            self.gscale = arrset(self.gscale,ig,float(-1.0))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),0.0)
        self.xupper = np.full((self.n,1),1.0)
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(0.2))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftp = []
        elftp = loaset(elftp,it,0,'PIJ')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['1']),int(v_['L'])+1):
            for J in range(int(v_['1']),int(v_['NN'])+1):
                ename = 'E'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eSQ')
                ielftype = arrset(ielftype,ie,iet_["eSQ"])
                vname = 'X'+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.0),float(1.0),float(0.2))
                posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                posep = np.where(elftp[ielftype[ie]]=='PIJ')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['P'+str(I)+','+str(J)]))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gNEXP',igt_)
        [it,igt_,_] = s2mpj_ii('gNEXP',igt_)
        grftp = []
        grftp = loaset(grftp,it,0,'CI')
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        self.grpar   = []
        for I in range(int(v_['1']),int(v_['L'])+1):
            ig = ig_['OBJ'+str(I)]
            self.grftype = arrset(self.grftype,ig,'gNEXP')
            for J in range(int(v_['1']),int(v_['NN'])+1):
                ig = ig_['OBJ'+str(I)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(I)+','+str(J)])
                self.grelw = loaset(self.grelw,ig,posel,float(v_['A'+str(I)+','+str(J)]))
            ig = ig_['OBJ'+str(I)]
            posgp = np.where(grftp[igt_[self.grftype[ig]]]=='CI')[0]
            self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['C'+str(I)]))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               -3.32288689158
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass   = "C-COBR2-AN-6-0"
        self.objderlvl = 2

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSQ(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = (EV_[0]-self.elpar[iel_][0])*(EV_[0]-self.elpar[iel_][0])
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*(EV_[0]-self.elpar[iel_][0])
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
    def gNEXP(self,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= self.grpar[igr_][0]*np.exp(-GVAR_)
        if nargout>1:
            g_ = -self.grpar[igr_][0]*np.exp(-GVAR_)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = self.grpar[igr_][0]*np.exp(-GVAR_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

