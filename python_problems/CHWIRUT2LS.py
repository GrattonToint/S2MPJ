from s2mpjlib import *
class  CHWIRUT2LS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CHWIRUT2LS
#    *********
# 
#    NIST Data fitting problem CHWIRUT2.
# 
#    Fit: y = exp[-b1*x]/(b2+b3*x) + e
# 
#    Source:  Problem from the NIST nonlinear regression test set
#      http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
# 
#    Reference: Chwirut, D., NIST (197?).  
#      Ultrasonic Reference Block Study. 
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

    name = 'CHWIRUT2LS'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['M'] = 54
        v_['N'] = 3
        v_['1'] = 1
        v_['X1'] = 0.5
        v_['X2'] = 1.0
        v_['X3'] = 1.75
        v_['X4'] = 3.75
        v_['X5'] = 5.75
        v_['X6'] = 0.875
        v_['X7'] = 2.25
        v_['X8'] = 3.25
        v_['X9'] = 5.25
        v_['X10'] = 0.75
        v_['X11'] = 1.75
        v_['X12'] = 2.75
        v_['X13'] = 4.75
        v_['X14'] = 0.625
        v_['X15'] = 1.25
        v_['X16'] = 2.25
        v_['X17'] = 4.25
        v_['X18'] = 0.5
        v_['X19'] = 3.0
        v_['X20'] = 0.75
        v_['X21'] = 3.0
        v_['X22'] = 1.5
        v_['X23'] = 6.0
        v_['X24'] = 3.0
        v_['X25'] = 6.0
        v_['X26'] = 1.5
        v_['X27'] = 3.0
        v_['X28'] = 0.5
        v_['X29'] = 2.0
        v_['X30'] = 4.0
        v_['X31'] = 0.75
        v_['X32'] = 2.0
        v_['X33'] = 5.0
        v_['X34'] = 0.75
        v_['X35'] = 2.25
        v_['X36'] = 3.75
        v_['X37'] = 5.75
        v_['X38'] = 3.0
        v_['X39'] = 0.75
        v_['X40'] = 2.5
        v_['X41'] = 4.0
        v_['X42'] = 0.75
        v_['X43'] = 2.5
        v_['X44'] = 4.0
        v_['X45'] = 0.75
        v_['X46'] = 2.5
        v_['X47'] = 4.0
        v_['X48'] = 0.5
        v_['X49'] = 6.0
        v_['X50'] = 3.0
        v_['X51'] = 0.5
        v_['X52'] = 2.75
        v_['X53'] = 0.5
        v_['X54'] = 1.75
        v_['Y1'] = 92.9
        v_['Y2'] = 57.1
        v_['Y3'] = 31.05
        v_['Y4'] = 11.5875
        v_['Y5'] = 8.025
        v_['Y6'] = 63.6
        v_['Y7'] = 21.4
        v_['Y8'] = 14.25
        v_['Y9'] = 8.475
        v_['Y10'] = 63.8
        v_['Y11'] = 26.8
        v_['Y12'] = 16.4625
        v_['Y13'] = 7.125
        v_['Y14'] = 67.3
        v_['Y15'] = 41.0
        v_['Y16'] = 21.15
        v_['Y17'] = 8.175
        v_['Y18'] = 81.50
        v_['Y19'] = 13.12
        v_['Y20'] = 59.9
        v_['Y21'] = 14.62
        v_['Y22'] = 32.9
        v_['Y23'] = 5.44
        v_['Y24'] = 12.56
        v_['Y25'] = 5.44
        v_['Y26'] = 32.0
        v_['Y27'] = 13.95
        v_['Y28'] = 75.8
        v_['Y29'] = 20.0
        v_['Y30'] = 10.42
        v_['Y31'] = 59.5
        v_['Y32'] = 21.67
        v_['Y33'] = 8.55
        v_['Y34'] = 62.0
        v_['Y35'] = 20.2
        v_['Y36'] = 7.76
        v_['Y37'] = 3.75
        v_['Y38'] = 11.81
        v_['Y39'] = 54.7
        v_['Y40'] = 23.7
        v_['Y41'] = 11.55
        v_['Y42'] = 61.3
        v_['Y43'] = 17.7
        v_['Y44'] = 8.74
        v_['Y45'] = 59.2
        v_['Y46'] = 16.3
        v_['Y47'] = 8.62
        v_['Y48'] = 81.0
        v_['Y49'] = 4.87
        v_['Y50'] = 14.62
        v_['Y51'] = 81.7
        v_['Y52'] = 17.17
        v_['Y53'] = 81.3
        v_['Y54'] = 28.9
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
        self.x0[ix_['B1']] = float(0.1)
        self.x0[ix_['B2']] = float(0.01)
        self.x0[ix_['B3']] = float(0.02)
        pass
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eE16', iet_)
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
            self.elftype = arrset(self.elftype,ie,'eE16')
            ielftype = arrset(ielftype,ie,iet_["eE16"])
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
    def eE16(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        E = np.exp(-EV_[0,0]*self.elpar[iel_][0])
        EX = E*self.elpar[iel_][0]
        EX2 = EX*self.elpar[iel_][0]
        V2PV3X = EV_[1,0]+EV_[2,0]*self.elpar[iel_][0]
        V2PV32 = V2PV3X*V2PV3X
        V2PV33 = V2PV3X*V2PV32
        f_   = E/V2PV3X
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -EX/V2PV3X
            g_[1] = -E/V2PV32
            g_[2] = -EX/V2PV32
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = EX2/V2PV3X
                H_[0,1] = EX/V2PV32
                H_[1,0] = H_[0,1]
                H_[0,2] = EX2/V2PV32
                H_[2,0] = H_[0,2]
                H_[1,1] = 2.0*E/V2PV33
                H_[1,2] = 2.0*EX/V2PV33
                H_[2,1] = H_[1,2]
                H_[2,2] = 2.0*EX2/V2PV33
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

