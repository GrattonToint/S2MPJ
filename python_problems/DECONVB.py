from s2mpjlib import *
class  DECONVB(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DECONVB
#    *********
# 
#    A problem arising in deconvolution analysis 
#    (bounded variables version).
# 
#    Source:  
#    J.P. Rasson, Private communication, 1996.
# 
#    SIF input: Ph. Toint, Nov 1996.
#    unititialized variables fixed at zero, Nick Gould, Feb, 2013
# 
#    classification = "C-CSBR2-MN-61-0"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'DECONVB'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['0'] = 0
        v_['1'] = 1
        v_['LGSG'] = 11
        v_['LGTR'] = 40
        v_['-LGSG'] = -1*v_['LGSG']
        v_['PIC'] = 3.0000000000
        v_['TR1'] = 0.0000000000
        v_['TR2'] = 0.0000000000
        v_['TR3'] = 1.600000E-03
        v_['TR4'] = 5.400000E-03
        v_['TR5'] = 7.020000E-02
        v_['TR6'] = 0.1876000000
        v_['TR7'] = 0.3320000000
        v_['TR8'] = 0.7640000000
        v_['TR9'] = 0.9320000000
        v_['TR10'] = 0.8120000000
        v_['TR11'] = 0.3464000000
        v_['TR12'] = 0.2064000000
        v_['TR13'] = 8.300000E-02
        v_['TR14'] = 3.400000E-02
        v_['TR15'] = 6.179999E-02
        v_['TR16'] = 1.2000000000
        v_['TR17'] = 1.8000000000
        v_['TR18'] = 2.4000000000
        v_['TR19'] = 9.0000000000
        v_['TR20'] = 2.4000000000
        v_['TR21'] = 1.8010000000
        v_['TR22'] = 1.3250000000
        v_['TR23'] = 7.620000E-02
        v_['TR24'] = 0.2104000000
        v_['TR25'] = 0.2680000000
        v_['TR26'] = 0.5520000000
        v_['TR27'] = 0.9960000000
        v_['TR28'] = 0.3600000000
        v_['TR29'] = 0.2400000000
        v_['TR30'] = 0.1510000000
        v_['TR31'] = 2.480000E-02
        v_['TR32'] = 0.2432000000
        v_['TR33'] = 0.3602000000
        v_['TR34'] = 0.4800000000
        v_['TR35'] = 1.8000000000
        v_['TR36'] = 0.4800000000
        v_['TR37'] = 0.3600000000
        v_['TR38'] = 0.2640000000
        v_['TR39'] = 6.000000E-03
        v_['TR40'] = 6.000000E-03
        v_['SSG1'] = 1.000000E-02
        v_['SSG2'] = 2.000000E-02
        v_['SSG3'] = 0.4000000000
        v_['SSG4'] = 0.6000000000
        v_['SSG5'] = 0.8000000000
        v_['SSG6'] = 3.0000000000
        v_['SSG7'] = 0.8000000000
        v_['SSG8'] = 0.6000000000
        v_['SSG9'] = 0.4400000000
        v_['SSG10'] = 1.000000E-02
        v_['SSG11'] = 1.000000E-02
        v_['CC1'] = 0.0
        v_['CC2'] = 0.0
        v_['CC3'] = 0.0
        v_['CC4'] = 0.0
        v_['CC5'] = 0.0
        v_['CC6'] = 0.0
        v_['CC7'] = 0.0
        v_['CC8'] = 0.0
        v_['CC9'] = 0.0
        v_['CC10'] = 0.0
        v_['CC11'] = 0.0
        v_['CC12'] = 0.0
        v_['CC13'] = 0.0
        v_['CC14'] = 0.0
        v_['CC15'] = 0.0
        v_['CC16'] = 0.0
        v_['CC17'] = 0.0
        v_['CC18'] = 0.0
        v_['CC19'] = 0.0
        v_['CC20'] = 0.0
        v_['CC21'] = 0.0
        v_['CC22'] = 0.0
        v_['CC23'] = 0.0
        v_['CC24'] = 0.0
        v_['CC25'] = 0.0
        v_['CC26'] = 0.0
        v_['CC27'] = 0.0
        v_['CC28'] = 0.0
        v_['CC29'] = 0.0
        v_['CC30'] = 0.0
        v_['CC31'] = 0.0
        v_['CC32'] = 0.0
        v_['CC33'] = 0.0
        v_['CC34'] = 0.0
        v_['CC35'] = 0.0
        v_['CC36'] = 0.0
        v_['CC37'] = 0.0
        v_['CC38'] = 0.0
        v_['CC39'] = 0.0
        v_['CC40'] = 0.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for K in range(int(v_['-LGSG']),int(v_['LGTR'])+1):
            [iv,ix_,_] = s2mpj_ii('C'+str(K),ix_)
            self.xnames=arrset(self.xnames,iv,'C'+str(K))
        for I in range(int(v_['1']),int(v_['LGSG'])+1):
            [iv,ix_,_] = s2mpj_ii('SG'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'SG'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for K in range(int(v_['1']),int(v_['LGTR'])+1):
            [ig,ig_,_] = s2mpj_ii('R'+str(K),ig_)
            gtype = arrset(gtype,ig,'<>')
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        selfnob      = ngrp
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        for K in range(int(v_['1']),int(v_['LGTR'])+1):
            self.gconst = arrset(self.gconst,ig_['R'+str(K)],float(v_['TR'+str(K)]))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        for I in range(int(v_['1']),int(v_['LGSG'])+1):
            self.xlower[ix_['SG'+str(I)]] = 0.0
            self.xupper[ix_['SG'+str(I)]] = v_['PIC']
        for K in range(int(v_['-LGSG']),int(v_['0'])+1):
            self.xlower[ix_['C'+str(K)]] = 0.0
            self.xupper[ix_['C'+str(K)]] = 0.0
        for K in range(int(v_['1']),int(v_['LGTR'])+1):
            self.xlower[ix_['C'+str(K)]] = 0.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        for K in range(int(v_['1']),int(v_['LGTR'])+1):
            self.x0[ix_['C'+str(K)]] = float(v_['CC'+str(K)])
        for I in range(int(v_['1']),int(v_['LGSG'])+1):
            self.x0[ix_['SG'+str(I)]] = float(v_['SSG'+str(I)])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'ePR', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        elftp = []
        elftp = loaset(elftp,it,0,'IDX')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for K in range(int(v_['1']),int(v_['LGTR'])+1):
            for I in range(int(v_['1']),int(v_['LGSG'])+1):
                v_['K-I'] = K-I
                v_['K-I+1'] = 1+v_['K-I']
                v_['RIDX'] = float(v_['K-I+1'])
                ename = 'PROD'+str(K)+','+str(I)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'ePR')
                ielftype = arrset(ielftype,ie,iet_["ePR"])
                vname = 'SG'+str(I)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                posev = np.where(elftv[ielftype[ie]]=='X')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'C'+str(int(v_['K-I+1']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                posev = np.where(elftv[ielftype[ie]]=='Y')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                posep = np.where(elftp[ielftype[ie]]=='IDX')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['RIDX']))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gSQ',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for K in range(int(v_['1']),int(v_['LGTR'])+1):
            ig = ig_['R'+str(K)]
            self.grftype = arrset(self.grftype,ig,'gSQ')
            for I in range(int(v_['1']),int(v_['LGSG'])+1):
                ig = ig_['R'+str(K)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['PROD'+str(K)+','+str(I)])
                self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass   = "C-CSBR2-MN-61-0"
        self.objderlvl = 2


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def ePR(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        NEGIDX = self.elpar[iel_][0]<=0.0
        if NEGIDX!=0:
            SCAL = 0.0
        if NEGIDX==0:
            SCAL = 1.0
        f_   = SCAL*EV_[0,0]*EV_[1,0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = SCAL*EV_[1,0]
            g_[1] = SCAL*EV_[0,0]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = SCAL
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gSQ(self,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= GVAR_**2
        if nargout>1:
            g_ = 2*GVAR_
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 2
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

