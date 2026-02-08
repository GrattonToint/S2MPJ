from s2mpjlib import *
class  ECKERLE4LS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : ECKERLE4LS
#    *********
# 
#    NIST Data fitting problem ECKERLE4.
# 
#    Fit: y = (b1/b2) * exp[-0.5*((x-b3)/b2)**2] + e
# 
#    Source:  Problem from the NIST nonlinear regression test set
#      http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
# 
#    Reference: Eckerle, K., NIST (197?).  
#      Circular Interference Transmittance Study.
# 
#    SIF input: Nick Gould and Tyrone Rees, Oct 2015
# 
#    classification = "C-CSUR2-MN-3-0"
# 
#    Number of data values
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'ECKERLE4LS'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['M'] = 35
        v_['N'] = 3
        v_['1'] = 1
        v_['X1'] = 400.000000
        v_['X2'] = 405.000000
        v_['X3'] = 410.000000
        v_['X4'] = 415.000000
        v_['X5'] = 420.000000
        v_['X6'] = 425.000000
        v_['X7'] = 430.000000
        v_['X8'] = 435.000000
        v_['X9'] = 436.500000
        v_['X10'] = 438.000000
        v_['X11'] = 439.500000
        v_['X12'] = 441.000000
        v_['X13'] = 442.500000
        v_['X14'] = 444.000000
        v_['X15'] = 445.500000
        v_['X16'] = 447.000000
        v_['X17'] = 448.500000
        v_['X18'] = 450.000000
        v_['X19'] = 451.500000
        v_['X20'] = 453.000000
        v_['X21'] = 454.500000
        v_['X22'] = 456.000000
        v_['X23'] = 457.500000
        v_['X24'] = 459.000000
        v_['X25'] = 460.500000
        v_['X26'] = 462.000000
        v_['X27'] = 463.500000
        v_['X28'] = 465.000000
        v_['X29'] = 470.000000
        v_['X30'] = 475.000000
        v_['X31'] = 480.000000
        v_['X32'] = 485.000000
        v_['X33'] = 490.000000
        v_['X34'] = 495.000000
        v_['X35'] = 500.000000
        v_['Y1'] = 0.0001575
        v_['Y2'] = 0.0001699
        v_['Y3'] = 0.0002350
        v_['Y4'] = 0.0003102
        v_['Y5'] = 0.0004917
        v_['Y6'] = 0.0008710
        v_['Y7'] = 0.0017418
        v_['Y8'] = 0.0046400
        v_['Y9'] = 0.0065895
        v_['Y10'] = 0.0097302
        v_['Y11'] = 0.0149002
        v_['Y12'] = 0.0237310
        v_['Y13'] = 0.0401683
        v_['Y14'] = 0.0712559
        v_['Y15'] = 0.1264458
        v_['Y16'] = 0.2073413
        v_['Y17'] = 0.2902366
        v_['Y18'] = 0.3445623
        v_['Y19'] = 0.3698049
        v_['Y20'] = 0.3668534
        v_['Y21'] = 0.3106727
        v_['Y22'] = 0.2078154
        v_['Y23'] = 0.1164354
        v_['Y24'] = 0.0616764
        v_['Y25'] = 0.0337200
        v_['Y26'] = 0.0194023
        v_['Y27'] = 0.0117831
        v_['Y28'] = 0.0074357
        v_['Y29'] = 0.0022732
        v_['Y30'] = 0.0008800
        v_['Y31'] = 0.0004579
        v_['Y32'] = 0.0002345
        v_['Y33'] = 0.0001586
        v_['Y34'] = 0.0001143
        v_['Y35'] = 0.0000710
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
        self.x0[ix_['B1']] = float(1.0)
        self.x0[ix_['B2']] = float(10.0)
        self.x0[ix_['B3']] = float(500.0)
        pass
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eE', iet_)
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
            self.elftype = arrset(self.elftype,ie,'eE')
            ielftype = arrset(ielftype,ie,iet_["eE"])
            vname = 'B1'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'B2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'B3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
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
    def eE(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        V3MX = EV_[2,0]-self.elpar[iel_][0]
        TV3MX = 2.0*V3MX
        V3MX2 = V3MX**2
        V22 = EV_[1,0]**2
        V23 = EV_[1,0]*V22
        V24 = V22*V22
        V25 = V22*V23
        V26 = V23*V23
        V27 = V23*V24
        E = np.exp(-0.5*V3MX2/V22)
        V1E = EV_[0,0]*E
        DIFF = V3MX2/V24-1.0/V22
        f_   = EV_[0,0]*E/EV_[1,0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = E/EV_[1,0]
            g_[1] = V1E*DIFF
            g_[2] = -V1E*V3MX/V23
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = E*DIFF
                H_[1,0] = H_[0,1]
                H_[0,2] = -0.5*E*TV3MX/V23
                H_[2,0] = H_[0,2]
                H_[1,1] = 2.0*V1E/V23-5.0*V1E*V3MX2/V25+V1E*V3MX**4/V27
                H_[1,2] = 1.5*V1E*TV3MX/V24-0.5*V1E*TV3MX*V3MX2/V26
                H_[2,1] = H_[1,2]
                H_[2,2] = 0.5*V1E*V3MX2/V25-V1E/V23
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

