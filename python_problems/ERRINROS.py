from s2mpjlib import *
class  ERRINROS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : ERRINROS
#    --------
# 
#    A nonlinear function similar to the chained Rosenbrock
#    problem CHNROSNB.
# 
#    Source:
#    An error in specifying problem CHNROSNB.
#    SIF input: Ph. Toint, Sept 1990.
# 
#    classification = "C-CSUR2-AN-V-0"
# 
#    Number of variables (at most 50)
# 
#           Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER
# IE N                   25             $-PARAMETER
# IE N                   50             $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 21 VI 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'ERRINROS'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(10);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
        v_['ALPH1'] = 1.25
        v_['ALPH2'] = 1.40
        v_['ALPH3'] = 2.40
        v_['ALPH4'] = 1.40
        v_['ALPH5'] = 1.75
        v_['ALPH6'] = 1.20
        v_['ALPH7'] = 2.25
        v_['ALPH8'] = 1.20
        v_['ALPH9'] = 1.00
        v_['ALPH10'] = 1.10
        v_['ALPH11'] = 1.50
        v_['ALPH12'] = 1.60
        v_['ALPH13'] = 1.25
        v_['ALPH14'] = 1.25
        v_['ALPH15'] = 1.20
        v_['ALPH16'] = 1.20
        v_['ALPH17'] = 1.40
        v_['ALPH18'] = 0.50
        v_['ALPH19'] = 0.50
        v_['ALPH20'] = 1.25
        v_['ALPH21'] = 1.80
        v_['ALPH22'] = 0.75
        v_['ALPH23'] = 1.25
        v_['ALPH24'] = 1.40
        v_['ALPH25'] = 1.60
        v_['ALPH26'] = 2.00
        v_['ALPH27'] = 1.00
        v_['ALPH28'] = 1.60
        v_['ALPH29'] = 1.25
        v_['ALPH30'] = 2.75
        v_['ALPH31'] = 1.25
        v_['ALPH32'] = 1.25
        v_['ALPH33'] = 1.25
        v_['ALPH34'] = 3.00
        v_['ALPH35'] = 1.50
        v_['ALPH36'] = 2.00
        v_['ALPH37'] = 1.25
        v_['ALPH38'] = 1.40
        v_['ALPH39'] = 1.80
        v_['ALPH40'] = 1.50
        v_['ALPH41'] = 2.20
        v_['ALPH42'] = 1.40
        v_['ALPH43'] = 1.50
        v_['ALPH44'] = 1.25
        v_['ALPH45'] = 2.00
        v_['ALPH46'] = 1.50
        v_['ALPH47'] = 1.25
        v_['ALPH48'] = 1.40
        v_['ALPH49'] = 0.60
        v_['ALPH50'] = 1.50
        v_['1'] = 1
        v_['2'] = 2
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
        for I in range(int(v_['2']),int(v_['N'])+1):
            v_['I-1'] = -1+I
            [ig,ig_,_] = s2mpj_ii('SQ'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['I-1']))]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('B'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(I)]])
            valA = np.append(valA,float(1.0))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        for I in range(int(v_['2']),int(v_['N'])+1):
            self.gconst = arrset(self.gconst,ig_['B'+str(I)],float(1.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        for I in range(int(v_['1']),int(v_['N'])+1):
            self.x0[ix_['X'+str(I)]] = float(-1.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eETYPE', iet_)
        elftv = loaset(elftv,it,0,'V1')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['2']),int(v_['N'])+1):
            ename = 'ELA'+str(I)
            [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
            if newelt:
                self.elftype = arrset(self.elftype,ie,'eETYPE')
                ielftype = arrset(ielftype,ie,iet_['eETYPE'])
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
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
        for I in range(int(v_['2']),int(v_['N'])+1):
            v_['AI2'] = v_['ALPH'+str(I)]*v_['ALPH'+str(I)]
            v_['AI'] = 16.0*v_['AI2']
            ig = ig_['SQ'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['ELA'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,float(v_['AI']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 0.0
#    Solution
# LO SOLTN(10)            6.69463214
# LO SOLTN(25)            18.4609060
# LO SOLTN(50)            39.9041540
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass   = "C-CSUR2-AN-V-0"
        self.objderlvl = 2

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eETYPE(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = -EV_[0]**2
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -2.0*EV_[0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = -2.0
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

