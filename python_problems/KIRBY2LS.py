from s2mpjlib import *
class  KIRBY2LS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : KIRBY2LS
#    *********
# 
#    NIST Data fitting problem KIRBY2.
# 
#    Fit: y = (b1 + b2*x + b3*x**2) /(1 + b4*x + b5*x**2) + e
# 
#    Source:  Problem from the NIST nonlinear regression test set
#      http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
# 
#    Reference: Kirby, R., NIST (197?).  
#      Scanning electron microscope line width standards.
# 
#    SIF input: Nick Gould and Tyrone Rees, Oct 2015
# 
#    classification = "C-CSUR2-MN-5-0"
# 
#    Number of data values
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'KIRBY2LS'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['M'] = 151
        v_['N'] = 5
        v_['1'] = 1
        v_['X1'] = 9.65
        v_['X2'] = 10.74
        v_['X3'] = 11.81
        v_['X4'] = 12.88
        v_['X5'] = 14.06
        v_['X6'] = 15.28
        v_['X7'] = 16.63
        v_['X8'] = 18.19
        v_['X9'] = 19.88
        v_['X10'] = 21.84
        v_['X11'] = 24.0
        v_['X12'] = 26.25
        v_['X13'] = 28.86
        v_['X14'] = 31.85
        v_['X15'] = 35.79
        v_['X16'] = 40.18
        v_['X17'] = 44.74
        v_['X18'] = 49.53
        v_['X19'] = 53.94
        v_['X20'] = 58.29
        v_['X21'] = 62.63
        v_['X22'] = 67.03
        v_['X23'] = 71.25
        v_['X24'] = 75.22
        v_['X25'] = 79.33
        v_['X26'] = 83.56
        v_['X27'] = 87.75
        v_['X28'] = 91.93
        v_['X29'] = 96.1
        v_['X30'] = 100.28
        v_['X31'] = 104.46
        v_['X32'] = 108.66
        v_['X33'] = 112.71
        v_['X34'] = 116.88
        v_['X35'] = 121.33
        v_['X36'] = 125.79
        v_['X37'] = 125.79
        v_['X38'] = 128.74
        v_['X39'] = 130.27
        v_['X40'] = 133.33
        v_['X41'] = 134.79
        v_['X42'] = 137.93
        v_['X43'] = 139.33
        v_['X44'] = 142.46
        v_['X45'] = 143.9
        v_['X46'] = 146.91
        v_['X47'] = 148.51
        v_['X48'] = 151.41
        v_['X49'] = 153.17
        v_['X50'] = 155.97
        v_['X51'] = 157.76
        v_['X52'] = 160.56
        v_['X53'] = 162.30
        v_['X54'] = 165.21
        v_['X55'] = 166.9
        v_['X56'] = 169.92
        v_['X57'] = 170.32
        v_['X58'] = 171.54
        v_['X59'] = 173.79
        v_['X60'] = 174.57
        v_['X61'] = 176.25
        v_['X62'] = 177.34
        v_['X63'] = 179.19
        v_['X64'] = 181.02
        v_['X65'] = 182.08
        v_['X66'] = 183.88
        v_['X67'] = 185.75
        v_['X68'] = 186.80
        v_['X69'] = 188.63
        v_['X70'] = 190.45
        v_['X71'] = 191.48
        v_['X72'] = 193.35
        v_['X73'] = 195.22
        v_['X74'] = 196.23
        v_['X75'] = 198.05
        v_['X76'] = 199.97
        v_['X77'] = 201.06
        v_['X78'] = 202.83
        v_['X79'] = 204.69
        v_['X80'] = 205.86
        v_['X81'] = 207.58
        v_['X82'] = 209.50
        v_['X83'] = 210.65
        v_['X84'] = 212.33
        v_['X85'] = 215.43
        v_['X86'] = 217.16
        v_['X87'] = 220.21
        v_['X88'] = 221.98
        v_['X89'] = 225.06
        v_['X90'] = 226.79
        v_['X91'] = 229.92
        v_['X92'] = 231.69
        v_['X93'] = 234.77
        v_['X94'] = 236.6
        v_['X95'] = 239.63
        v_['X96'] = 241.50
        v_['X97'] = 244.48
        v_['X98'] = 246.40
        v_['X99'] = 249.35
        v_['X100'] = 251.32
        v_['X101'] = 254.22
        v_['X102'] = 256.24
        v_['X103'] = 259.11
        v_['X104'] = 261.18
        v_['X105'] = 264.02
        v_['X106'] = 266.13
        v_['X107'] = 268.94
        v_['X108'] = 271.09
        v_['X109'] = 273.87
        v_['X110'] = 276.08
        v_['X111'] = 278.83
        v_['X112'] = 281.08
        v_['X113'] = 283.81
        v_['X114'] = 286.11
        v_['X115'] = 288.81
        v_['X116'] = 291.08
        v_['X117'] = 293.75
        v_['X118'] = 295.99
        v_['X119'] = 298.64
        v_['X120'] = 300.84
        v_['X121'] = 302.02
        v_['X122'] = 303.48
        v_['X123'] = 305.65
        v_['X124'] = 308.27
        v_['X125'] = 310.41
        v_['X126'] = 313.01
        v_['X127'] = 315.12
        v_['X128'] = 317.71
        v_['X129'] = 319.79
        v_['X130'] = 322.36
        v_['X131'] = 324.42
        v_['X132'] = 326.98
        v_['X133'] = 329.01
        v_['X134'] = 331.56
        v_['X135'] = 333.56
        v_['X136'] = 336.1
        v_['X137'] = 338.08
        v_['X138'] = 340.6
        v_['X139'] = 342.57
        v_['X140'] = 345.08
        v_['X141'] = 347.02
        v_['X142'] = 349.52
        v_['X143'] = 351.44
        v_['X144'] = 353.93
        v_['X145'] = 355.83
        v_['X146'] = 358.32
        v_['X147'] = 360.2
        v_['X148'] = 362.67
        v_['X149'] = 364.53
        v_['X150'] = 367.0
        v_['X151'] = 371.3
        v_['Y1'] = 0.0082
        v_['Y2'] = 0.0112
        v_['Y3'] = 0.0149
        v_['Y4'] = 0.0198
        v_['Y5'] = 0.0248
        v_['Y6'] = 0.0324
        v_['Y7'] = 0.042
        v_['Y8'] = 0.0549
        v_['Y9'] = 0.0719
        v_['Y10'] = 0.0963
        v_['Y11'] = 0.1291
        v_['Y12'] = 0.171
        v_['Y13'] = 0.2314
        v_['Y14'] = 0.3227
        v_['Y15'] = 0.4809
        v_['Y16'] = 0.7084
        v_['Y17'] = 1.022
        v_['Y18'] = 1.458
        v_['Y19'] = 1.952
        v_['Y20'] = 2.541
        v_['Y21'] = 3.223
        v_['Y22'] = 3.999
        v_['Y23'] = 4.852
        v_['Y24'] = 5.732
        v_['Y25'] = 6.727
        v_['Y26'] = 7.835
        v_['Y27'] = 9.025
        v_['Y28'] = 10.267
        v_['Y29'] = 11.578
        v_['Y30'] = 12.944
        v_['Y31'] = 14.377
        v_['Y32'] = 15.856
        v_['Y33'] = 17.331
        v_['Y34'] = 18.885
        v_['Y35'] = 20.575
        v_['Y36'] = 22.32
        v_['Y37'] = 22.303
        v_['Y38'] = 23.46
        v_['Y39'] = 24.06
        v_['Y40'] = 25.272
        v_['Y41'] = 25.853
        v_['Y42'] = 27.11
        v_['Y43'] = 27.658
        v_['Y44'] = 28.924
        v_['Y45'] = 29.511
        v_['Y46'] = 30.71
        v_['Y47'] = 31.35
        v_['Y48'] = 32.52
        v_['Y49'] = 33.23
        v_['Y50'] = 34.33
        v_['Y51'] = 35.06
        v_['Y52'] = 36.17
        v_['Y53'] = 36.84
        v_['Y54'] = 38.01
        v_['Y55'] = 38.67
        v_['Y56'] = 39.87
        v_['Y57'] = 40.03
        v_['Y58'] = 40.5
        v_['Y59'] = 41.37
        v_['Y60'] = 41.67
        v_['Y61'] = 42.31
        v_['Y62'] = 42.73
        v_['Y63'] = 43.46
        v_['Y64'] = 44.14
        v_['Y65'] = 44.55
        v_['Y66'] = 45.22
        v_['Y67'] = 45.92
        v_['Y68'] = 46.3
        v_['Y69'] = 47.0
        v_['Y70'] = 47.68
        v_['Y71'] = 48.06
        v_['Y72'] = 48.74
        v_['Y73'] = 49.41
        v_['Y74'] = 49.76
        v_['Y75'] = 50.43
        v_['Y76'] = 51.11
        v_['Y77'] = 51.5
        v_['Y78'] = 52.12
        v_['Y79'] = 52.76
        v_['Y80'] = 53.18
        v_['Y81'] = 53.78
        v_['Y82'] = 54.46
        v_['Y83'] = 54.83
        v_['Y84'] = 55.4
        v_['Y85'] = 56.43
        v_['Y86'] = 57.03
        v_['Y87'] = 58.0
        v_['Y88'] = 58.61
        v_['Y89'] = 59.58
        v_['Y90'] = 60.11
        v_['Y91'] = 61.1
        v_['Y92'] = 61.65
        v_['Y93'] = 62.59
        v_['Y94'] = 63.12
        v_['Y95'] = 64.03
        v_['Y96'] = 64.62
        v_['Y97'] = 65.49
        v_['Y98'] = 66.03
        v_['Y99'] = 66.89
        v_['Y100'] = 67.42
        v_['Y101'] = 68.23
        v_['Y102'] = 68.77
        v_['Y103'] = 69.59
        v_['Y104'] = 70.11
        v_['Y105'] = 70.86
        v_['Y106'] = 71.43
        v_['Y107'] = 72.16
        v_['Y108'] = 72.7
        v_['Y109'] = 73.4
        v_['Y110'] = 73.93
        v_['Y111'] = 74.6
        v_['Y112'] = 75.16
        v_['Y113'] = 75.82
        v_['Y114'] = 76.34
        v_['Y115'] = 76.98
        v_['Y116'] = 77.48
        v_['Y117'] = 78.08
        v_['Y118'] = 78.6
        v_['Y119'] = 79.17
        v_['Y120'] = 79.62
        v_['Y121'] = 79.88
        v_['Y122'] = 80.19
        v_['Y123'] = 80.66
        v_['Y124'] = 81.22
        v_['Y125'] = 81.66
        v_['Y126'] = 82.16
        v_['Y127'] = 82.59
        v_['Y128'] = 83.14
        v_['Y129'] = 83.5
        v_['Y130'] = 84.0
        v_['Y131'] = 84.4
        v_['Y132'] = 84.89
        v_['Y133'] = 85.26
        v_['Y134'] = 85.74
        v_['Y135'] = 86.07
        v_['Y136'] = 86.54
        v_['Y137'] = 86.89
        v_['Y138'] = 87.32
        v_['Y139'] = 87.65
        v_['Y140'] = 88.1
        v_['Y141'] = 88.43
        v_['Y142'] = 88.83
        v_['Y143'] = 89.12
        v_['Y144'] = 89.54
        v_['Y145'] = 89.85
        v_['Y146'] = 90.25
        v_['Y147'] = 90.55
        v_['Y148'] = 90.93
        v_['Y149'] = 91.2
        v_['Y150'] = 91.55
        v_['Y151'] = 92.2
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
        self.x0[ix_['B1']] = float(2.0)
        self.x0[ix_['B2']] = float(-0.1)
        self.x0[ix_['B3']] = float(0.003)
        self.x0[ix_['B4']] = float(-0.001)
        self.x0[ix_['B5']] = float(0.00001)
        pass
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eE18', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftv = loaset(elftv,it,3,'V4')
        elftv = loaset(elftv,it,4,'V5')
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
            self.elftype = arrset(self.elftype,ie,'eE18')
            ielftype = arrset(ielftype,ie,iet_["eE18"])
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
            vname = 'B4'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V4')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'B5'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V5')[0]
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
        self.pbclass   = "C-CSUR2-MN-5-0"
        self.objderlvl = 2

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eE18(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        X2 = self.elpar[iel_][0]*self.elpar[iel_][0]
        X3 = self.elpar[iel_][0]*X2
        X4 = self.elpar[iel_][0]*X3
        T = EV_[0]+EV_[1]*self.elpar[iel_][0]+EV_[2]*X2
        B = 1.0+EV_[3]*self.elpar[iel_][0]+EV_[4]*X2
        B2 = B*B
        B3 = B*B2
        f_   = T/B
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 1/B
            g_[1] = self.elpar[iel_][0]/B
            g_[2] = X2/B
            g_[3] = -T*self.elpar[iel_][0]/B2
            g_[4] = -T*X2/B2
            if nargout>2:
                H_ = np.zeros((5,5))
                H_[0,3] = -self.elpar[iel_][0]/B2
                H_[3,0] = H_[0,3]
                H_[0,4] = -X2/B2
                H_[4,0] = H_[0,4]
                H_[1,3] = -X2/B2
                H_[3,1] = H_[1,3]
                H_[1,4] = -X3/B2
                H_[4,1] = H_[1,4]
                H_[2,3] = -X3/B2
                H_[3,2] = H_[2,3]
                H_[2,4] = -X4/B2
                H_[4,2] = H_[2,4]
                H_[3,3] = 2.0*T*X2/B3
                H_[3,4] = 2.0*T*X3/B3
                H_[4,3] = H_[3,4]
                H_[4,4] = 2.0*T*X4/B3
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

