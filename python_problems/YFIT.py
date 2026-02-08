from s2mpjlib import *
class  YFIT(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    A nonlinear least-squares problem.
# 
#    Source:
#    an exercize for L. Watson course on LANCELOT in the Spring 1993.
# 
#    SIF input: Brian E. Lindholm, Virginia Tech., Spring 1993.
#               derivatives corrected by Nick Gould, June 2019.
# 
#    classification = "C-CSBR2-MN-3-0"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'YFIT'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['zero'] = 0
        v_['p'] = 16
        v_['realp'] = 16.0
        v_['y0'] = 21.158931
        v_['y1'] = 17.591719
        v_['y2'] = 14.046854
        v_['y3'] = 10.519732
        v_['y4'] = 7.0058392
        v_['y5'] = 3.5007293
        v_['y6'] = 0.0000000
        v_['y7'] = -3.5007293
        v_['y8'] = -7.0058392
        v_['y9'] = -10.519732
        v_['y10'] = -14.046854
        v_['y11'] = -17.591719
        v_['y12'] = -21.158931
        v_['y13'] = -24.753206
        v_['y14'] = -28.379405
        v_['y15'] = -32.042552
        v_['y16'] = -35.747869
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [iv,ix_,_] = s2mpj_ii('alpha',ix_)
        self.xnames=arrset(self.xnames,iv,'alpha')
        [iv,ix_,_] = s2mpj_ii('beta',ix_)
        self.xnames=arrset(self.xnames,iv,'beta')
        [iv,ix_,_] = s2mpj_ii('dist',ix_)
        self.xnames=arrset(self.xnames,iv,'dist')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for i in range(int(v_['zero']),int(v_['p'])+1):
            [ig,ig_,_] = s2mpj_ii('diff'+str(i),ig_)
            gtype = arrset(gtype,ig,'<>')
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        selfnob      = ngrp
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        for i in range(int(v_['zero']),int(v_['p'])+1):
            self.gconst = arrset(self.gconst,ig_['diff'+str(i)],float(v_['y'+str(i)]))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        self.xlower[ix_['alpha']] = -float('Inf')
        self.xupper[ix_['alpha']] = +float('Inf')
        self.xlower[ix_['beta']] = -float('Inf')
        self.xupper[ix_['beta']] = +float('Inf')
        self.xlower[ix_['dist']] = 0.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.x0[ix_['alpha']] = float(0.60)
        self.x0[ix_['beta']] = float(-0.60)
        self.x0[ix_['dist']] = float(20.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'etanab', iet_)
        elftv = loaset(elftv,it,0,'a1')
        elftv = loaset(elftv,it,1,'b1')
        elftv = loaset(elftv,it,2,'d1')
        elftp = []
        elftp = loaset(elftp,it,0,'point')
        elftp = loaset(elftp,it,1,'count')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for i in range(int(v_['zero']),int(v_['p'])+1):
            v_['index'] = float(i)
            ename = 'est'+str(i)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'etanab')
            ielftype = arrset(ielftype,ie,iet_["etanab"])
            vname = 'alpha'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='a1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'beta'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='b1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'dist'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='d1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='point')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['index']))
            posep = np.where(elftp[ielftype[ie]]=='count')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['realp']))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gsquare',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for i in range(int(v_['zero']),int(v_['p'])+1):
            ig = ig_['diff'+str(i)]
            self.grftype = arrset(self.grftype,ig,'gsquare')
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['est'+str(i)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
# LO SOLUTION            0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass   = "C-CSBR2-MN-3-0"
        self.objderlvl = 2

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def etanab(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        frac = self.elpar[iel_][0]/self.elpar[iel_][1]
        ttan = np.tan(EV_[0,0]*(1.0-frac)+EV_[1,0]*frac)
        tsec = 1.0/np.cos(EV_[0,0]*(1.0-frac)+EV_[1,0]*frac)
        tsec2 = tsec*tsec
        f_   = EV_[2,0]*ttan
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[2,0]*(1.0-frac)*tsec2
            g_[1] = EV_[2,0]*frac*tsec2
            g_[2] = ttan
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = 2.0*EV_[2,0]*((1.0-frac)**2)*tsec2*ttan
                H_[1,1] = 2.0*EV_[2,0]*(frac**2)*tsec2*ttan
                H_[0,1] = 2.0*EV_[2,0]*(1.0-frac)*frac*tsec2*ttan
                H_[1,0] = H_[0,1]
                H_[0,2] = (1.0-frac)*tsec2
                H_[2,0] = H_[0,2]
                H_[1,2] = frac*tsec2
                H_[2,1] = H_[1,2]
                H_[2,2] = 0.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gsquare(self,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= GVAR_*GVAR_
        if nargout>1:
            g_ = 2.0*GVAR_
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

