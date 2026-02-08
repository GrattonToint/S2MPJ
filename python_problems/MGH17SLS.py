from s2mpjlib import *
class  MGH17SLS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : MGH17SLS
#    *********
# 
#    NIST Data fitting problem MGH17 given as an inconsistent set of
#    nonlinear equations (scaled version).
# 
#    Fit: y = b1 + b2*exp[-x*0.01*b4] + b3*exp[-x*0.01*b5] + e
# 
#    Source:  Problem from the NIST nonlinear regression test set
#      http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
# 
#    Reference: Osborne, M. R. (1972).
#     Some aspects of nonlinear least squares calculations.
#     In Numerical Methods for Nonlinear Optimization, Lootsma (Ed).
#     New York, NY:  Academic Press, pp. 171-189.
# 
#    SIF input: Nick Gould and Tyrone Rees, Oct 2015
#    Least-squares version of MGH17S.SIF, Nick Gould, Jan 2020
# 
#    classification = "C-CSUR2-MN-5-33"
# 
#    Number of data values
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'MGH17SLS'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['M'] = 33
        v_['N'] = 5
        v_['1'] = 1
        v_['X1'] = 0.0E+0
        v_['X2'] = 1.0E+1
        v_['X3'] = 2.0E+1
        v_['X4'] = 3.0E+1
        v_['X5'] = 4.0E+1
        v_['X6'] = 5.0E+1
        v_['X7'] = 6.0E+1
        v_['X8'] = 7.0E+1
        v_['X9'] = 8.0E+1
        v_['X10'] = 9.0E+1
        v_['X11'] = 1.0E+2
        v_['X12'] = 1.1E+2
        v_['X13'] = 1.2E+2
        v_['X14'] = 1.3E+2
        v_['X15'] = 1.4E+2
        v_['X16'] = 1.5E+2
        v_['X17'] = 1.6E+2
        v_['X18'] = 1.7E+2
        v_['X19'] = 1.8E+2
        v_['X20'] = 1.9E+2
        v_['X21'] = 2.0E+2
        v_['X22'] = 2.1E+2
        v_['X23'] = 2.2E+2
        v_['X24'] = 2.3E+2
        v_['X25'] = 2.4E+2
        v_['X26'] = 2.5E+2
        v_['X27'] = 2.6E+2
        v_['X28'] = 2.7E+2
        v_['X29'] = 2.8E+2
        v_['X30'] = 2.9E+2
        v_['X31'] = 3.0E+2
        v_['X32'] = 3.1E+2
        v_['X33'] = 3.2E+2
        v_['Y1'] = 8.44E-1
        v_['Y2'] = 9.08E-1
        v_['Y3'] = 9.32E-1
        v_['Y4'] = 9.36E-1
        v_['Y5'] = 9.25E-1
        v_['Y6'] = 9.08E-1
        v_['Y7'] = 8.81E-1
        v_['Y8'] = 8.50E-1
        v_['Y9'] = 8.18E-1
        v_['Y10'] = 7.84E-1
        v_['Y11'] = 7.51E-1
        v_['Y12'] = 7.18E-1
        v_['Y13'] = 6.85E-1
        v_['Y14'] = 6.58E-1
        v_['Y15'] = 6.28E-1
        v_['Y16'] = 6.03E-1
        v_['Y17'] = 5.80E-1
        v_['Y18'] = 5.58E-1
        v_['Y19'] = 5.38E-1
        v_['Y20'] = 5.22E-1
        v_['Y21'] = 5.06E-1
        v_['Y22'] = 4.90E-1
        v_['Y23'] = 4.78E-1
        v_['Y24'] = 4.67E-1
        v_['Y25'] = 4.57E-1
        v_['Y26'] = 4.48E-1
        v_['Y27'] = 4.38E-1
        v_['Y28'] = 4.31E-1
        v_['Y29'] = 4.24E-1
        v_['Y30'] = 4.20E-1
        v_['Y31'] = 4.14E-1
        v_['Y32'] = 4.11E-1
        v_['Y33'] = 4.06E-1
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
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['B1']])
            valA = np.append(valA,float(1.0))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        selfnob      = ngrp
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        for I in range(int(v_['1']),int(v_['M'])+1):
            self.gconst = arrset(self.gconst,ig_['F'+str(I)],float(v_['Y'+str(I)]))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.x0[ix_['B1']] = float(50.0)
        self.x0[ix_['B2']] = float(150.0)
        self.x0[ix_['B3']] = float(-100.0)
        self.x0[ix_['B4']] = float(100.0)
        self.x0[ix_['B5']] = float(200.0)
        pass
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eE2', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftp = []
        elftp = loaset(elftp,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['1']),int(v_['M'])+1):
            ename = 'EA'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eE2')
            ielftype = arrset(ielftype,ie,iet_["eE2"])
            vname = 'B2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'B4'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='X')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['X'+str(I)]))
            ename = 'EB'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eE2')
            ielftype = arrset(ielftype,ie,iet_["eE2"])
            vname = 'B3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'B5'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
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
            self.grelt = loaset(self.grelt,ig,posel,ie_['EA'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['EB'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        self.objlower = 0.0
#    Solution
# LO SOLTN
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass   = "C-CSUR2-MN-5-33"
        self.objderlvl = 2

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eE2(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        S = 0.01
        SX = S*self.elpar[iel_][0]
        E = np.exp(-EV_[1,0]*SX)
        f_   = EV_[0,0]*E
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = E
            g_[1] = -EV_[0,0]*SX*E
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = -SX*E
                H_[1,0] = H_[0,1]
                H_[1,1] = EV_[0,0]*SX*SX*E
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

