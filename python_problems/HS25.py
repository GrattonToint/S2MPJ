from s2mpjlib import *
class  HS25(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS25
#    *********
# 
#    A nonlinear least squares problem with bounds.
# 
#    Source: problem 25 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: J-M Collin, Mar 1990.
# 
#    classification = "C-CSBR2-AN-3-0"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS25'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 3
        v_['1'] = 1
        v_['99'] = 99
        v_['2/3'] = 0.6666666666
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
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['1']),int(v_['99'])+1):
            [ig,ig_,_] = s2mpj_ii('O'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        selfnob      = ngrp
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        for I in range(int(v_['1']),int(v_['99'])+1):
            v_['IR'] = float(I)
            v_['I/100'] = 0.01*v_['IR']
            self.gconst = arrset(self.gconst,ig_['O'+str(I)],float(v_['I/100']))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        self.xlower[ix_['X1']] = 0.1
        self.xupper[ix_['X1']] = 100.0
        self.xlower[ix_['X2']] = 0.0
        self.xupper[ix_['X2']] = 25.6
        self.xlower[ix_['X3']] = 0.0
        self.xupper[ix_['X3']] = 5.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.x0[ix_['X1']] = float(100.0)
        self.x0[ix_['X2']] = float(12.5)
        self.x0[ix_['X3']] = float(3.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eWFI', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        elftv = loaset(elftv,it,2,'Z')
        elftp = []
        elftp = loaset(elftp,it,0,'W')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['1']),int(v_['99'])+1):
            v_['IR'] = float(I)
            v_['I/100'] = 0.01*v_['IR']
            v_['LOG01I'] = np.log(v_['I/100'])
            v_['M50LOG'] = -50.0*v_['LOG01I']
            v_['EXPLOG'] = np.log(v_['M50LOG'])
            v_['EXPL2/3'] = v_['EXPLOG']*v_['2/3']
            v_['EXP2/3'] = np.exp(v_['EXPL2/3'])
            v_['UI'] = 25.0+v_['EXP2/3']
            ename = 'E'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eWFI')
            ielftype = arrset(ielftype,ie,iet_["eWFI"])
            vname = 'X1'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='Y')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='Z')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='W')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['UI']))
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
        for I in range(int(v_['1']),int(v_['99'])+1):
            ig = ig_['O'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 0.0
#    Solution
# LO SOLTN               0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass   = "C-CSBR2-AN-3-0"
        self.objderlvl = 2

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eWFI(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        XI = 1.0/EV_[0,0]
        X2I = XI*XI
        X3I = X2I*XI
        WMY = self.elpar[iel_][0]-EV_[1,0]
        WMYEZ = WMY**EV_[2,0]
        LWMY = np.log(WMY)
        EXPO = np.exp(-XI*WMYEZ)
        f_   = EXPO
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = X2I*WMYEZ*EXPO
            g_[1] = XI*EV_[2,0]*WMY**(EV_[2,0]-1.0)*EXPO
            g_[2] = -XI*LWMY*WMYEZ*EXPO
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = EXPO*WMYEZ*X3I*(-2.0+XI*WMY**EV_[2,0])
                H_[0,1] = EXPO*EV_[2,0]*X2I*WMY**(EV_[2,0]-1.0)*(-1.0+XI*WMYEZ)
                H_[1,0] = H_[0,1]
                H_[0,2] = EXPO*X2I*WMYEZ*LWMY*(1.0-XI*WMYEZ)
                H_[2,0] = H_[0,2]
                H_[1,1]  = (
                      EXPO*XI*WMY**(EV_[2,0]-2.0)*EV_[2,0]*(-EV_[2,0]+1.0+XI*EV_[2,0]*WMYEZ))
                H_[1,2] = EXPO*XI*WMY**(EV_[2,0]-1.0)*(1.0+EV_[2,0]*LWMY*(1.0-XI*WMYEZ))
                H_[2,1] = H_[1,2]
                H_[2,2] = EXPO*WMYEZ*XI*LWMY**2*(-1.0+XI*WMYEZ)
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

