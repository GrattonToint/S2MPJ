from s2mpjlib import *
class  MGH10SLS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : MGH10SLS
#    *********
# 
#    NIST Data fitting problem MGH10 given as an inconsistent set of
#    nonlinear equations (scaled version)
# 
#    Fit: y = 0.01 * b1 * exp[ 1000 * b2 / (x + 100 * b3 ) ] + e
# 
#    Source:  Problem from the NIST nonlinear regression test set
#      http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
# 
#    Reference: Meyer, R. R. (1970).
#      Theoretical and computational aspects of nonlinear
#      regression.  In Nonlinear Programming, Rosen,
#      Mangasarian and Ritter (Eds).
#      New York, NY: Academic Press, pp. 465-486.
# 
#    SIF input: Nick Gould and Tyrone Rees, Oct 2015
# 
#    classification = "C-CSUR2-MN-3-0"
#    Least-squares version of MGH10S.SIF, Nick Gould, Jan 2020
# 
#    Number of data values
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 25 XI 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'MGH10SLS'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['M'] = 16
        v_['N'] = 3
        v_['1'] = 1
        v_['X1'] = 5.00E+01
        v_['X2'] = 5.50E+01
        v_['X3'] = 6.00E+01
        v_['X4'] = 6.50E+01
        v_['X5'] = 7.00E+01
        v_['X6'] = 7.50E+01
        v_['X7'] = 8.00E+01
        v_['X8'] = 8.50E+01
        v_['X9'] = 9.00E+01
        v_['X10'] = 9.50E+01
        v_['X11'] = 1.00E+02
        v_['X12'] = 1.05E+02
        v_['X13'] = 1.10E+02
        v_['X14'] = 1.15E+02
        v_['X15'] = 1.20E+02
        v_['X16'] = 1.25E+02
        v_['Y1'] = 3.478e+04
        v_['Y2'] = 2.861e+04
        v_['Y3'] = 2.365e+04
        v_['Y4'] = 1.963e+04
        v_['Y5'] = 1.637e+04
        v_['Y6'] = 1.372e+04
        v_['Y7'] = 1.154e+04
        v_['Y8'] = 9.744e+03
        v_['Y9'] = 8.261e+03
        v_['Y10'] = 7.030e+03
        v_['Y11'] = 6.005e+03
        v_['Y12'] = 5.147e+03
        v_['Y13'] = 4.427e+03
        v_['Y14'] = 3.820e+03
        v_['Y15'] = 3.307e+03
        v_['Y16'] = 2.872e+03
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('B'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'B'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['1']),int(v_['M'])+1):
            [ig,ig_,_] = s2mpj_ii('F'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        for I in range(int(v_['1']),int(v_['M'])+1):
            self.gconst = arrset(self.gconst,ig_['F'+str(I)],float(v_['Y'+str(I)]))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.x0[ix_['B1']] = float(200.0)
        self.x0[ix_['B2']] = float(400.0)
        self.x0[ix_['B3']] = float(250.0)
        pass
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eE12', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftp = []
        elftp = loaset(elftp,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['1']),int(v_['M'])+1):
            ename = 'E'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eE12')
            ielftype = arrset(ielftype,ie,iet_["eE12"])
            vname = 'B1'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'B2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'B3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='X')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['X'+str(I)]))
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
        for I in range(int(v_['1']),int(v_['M'])+1):
            ig = ig_['F'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        self.objlower = 0.0
#    Solution
# LO SOLTN
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass   = "C-CSUR2-MN-3-0"
        self.objderlvl = 2

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eE12(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        S1 = 0.01
        S2 = 1000.0
        S3 = 100.0
        XPV3 = self.elpar[iel_][0]+S3*EV_[2]
        XPV32 = XPV3*XPV3
        XPV33 = XPV3*XPV32
        E = S1*np.exp(S2*EV_[1]/XPV3)
        F = EV_[0]*E
        f_   = F
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = E
            g_[1] = S2*F/XPV3
            g_[2] = -S2*S3*F*EV_[1]/XPV32
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = S2*E/XPV3
                H_[1,0] = H_[0,1]
                H_[0,2] = -S2*S3*EV_[1]*E/XPV32
                H_[2,0] = H_[0,2]
                H_[1,1] = S2*S2*F/XPV32
                H_[1,2] = -S2*S3*F*(1.0/XPV32+S2*EV_[1]/XPV33)
                H_[2,1] = H_[1,2]
                H_[2,2] = S2*S3*S3*F*EV_[1]*(S2*EV_[1]/XPV3**4+2.0/XPV33)
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

