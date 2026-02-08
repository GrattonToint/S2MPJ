from s2mpjlib import *
class  CLNLBEAM(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CLNLBEAM
#    *********
# 
#    An optimal control version of the CLamped NonLinear BEAM problem.
#    The energy of a beam of length 1 compressed by a force P is to be
#    minimized.  The control variable is the derivative of the deflection angle.
# 
#    The problem is discretized using the trapezoidal rule. It is non-convex.
# 
#    Source:
#    H. Maurer and H.D. Mittelman,
#    "The non-linear beam via optimal control with bound state variables",
#    Optimal Control Applications and Methods 12, pp. 19-31, 1991.
# 
#    SIF input: Ph. Toint, Nov 1993.
# 
#    classification = "C-COOR2-MN-V-V"
# 
#    Discretization: specify the number of interior points + 1
# 
#           Alternative values for the SIF file parameters:
# IE NI                  10             $-PARAMETER n=33, m=20
# IE NI                  50             $-PARAMETER n=153, m=100
# IE NI                  100            $-PARAMETER n=303, m=200
# IE NI                  500            $-PARAMETER n=1503, m=1000
# IE NI                  1000           $-PARAMETER n=3003, m=2000 original value
# IE NI                  2000           $-PARAMETER n=6003, m=4000
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'CLNLBEAM'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['NI'] = int(10);  #  SIF file default value
        else:
            v_['NI'] = int(args[0])
# IE NI                  5000           $-PARAMETER n=15003, m=10000
        if nargin<2:
            v_['ALPHA'] = float(350.0);  #  SIF file default value
        else:
            v_['ALPHA'] = float(args[1])
        v_['RNI'] = float(v_['NI'])
        v_['NI-1'] = -1+v_['NI']
        v_['H'] = 1.0/v_['RNI']
        v_['H/4'] = 0.25*v_['H']
        v_['H/2'] = 0.5*v_['H']
        v_['AH'] = v_['ALPHA']*v_['H']
        v_['AH/2'] = 0.5*v_['AH']
        v_['-H/2'] = -0.5*v_['H']
        v_['0'] = 0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['0']),int(v_['NI'])+1):
            [iv,ix_,_] = s2mpj_ii('T'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'T'+str(I))
        for I in range(int(v_['0']),int(v_['NI'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I))
        for I in range(int(v_['0']),int(v_['NI'])+1):
            [iv,ix_,_] = s2mpj_ii('U'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'U'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('ENERGY',ig_)
        gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['0']),int(v_['NI-1'])+1):
            v_['I+1'] = 1+I
            [ig,ig_,_] = s2mpj_ii('EX'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'EX'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['I+1']))]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(I)]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('ET'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'ET'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['T'+str(int(v_['I+1']))]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['T'+str(I)]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['U'+str(int(v_['I+1']))]])
            valA = np.append(valA,float(v_['-H/2']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['U'+str(I)]])
            valA = np.append(valA,float(v_['-H/2']))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        legrps = np.where(gtype=='<=')[0]
        eqgrps = np.where(gtype=='==')[0]
        gegrps = np.where(gtype=='>=')[0]
        self.nle = len(legrps)
        self.neq = len(eqgrps)
        self.nge = len(gegrps)
        self.m   = self.nle+self.neq+self.nge
        self.congrps = np.concatenate((legrps,eqgrps,gegrps))
        self.cnames = cnames[self.congrps]
        self.nob = ngrp-self.m
        self.objgrps = np.where(gtype=='<>')[0]
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        for I in range(int(v_['0']),int(v_['NI'])+1):
            self.xlower[ix_['X'+str(I)]] = -0.05
            self.xupper[ix_['X'+str(I)]] = 0.05
        for I in range(int(v_['0']),int(v_['NI'])+1):
            self.xlower[ix_['T'+str(I)]] = -1.0
            self.xupper[ix_['T'+str(I)]] = 1.0
        self.xlower[ix_['X'+str(int(v_['0']))]] = 0.0
        self.xupper[ix_['X'+str(int(v_['0']))]] = 0.0
        self.xlower[ix_['X'+str(int(v_['NI']))]] = 0.0
        self.xupper[ix_['X'+str(int(v_['NI']))]] = 0.0
        self.xlower[ix_['T'+str(int(v_['0']))]] = 0.0
        self.xupper[ix_['T'+str(int(v_['0']))]] = 0.0
        self.xlower[ix_['T'+str(int(v_['NI']))]] = 0.0
        self.xupper[ix_['T'+str(int(v_['NI']))]] = 0.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        for I in range(int(v_['0']),int(v_['NI'])+1):
            v_['RI'] = float(I)
            v_['TT'] = v_['RI']*v_['H']
            v_['CTT'] = np.cos(v_['TT'])
            v_['SCTT'] = 0.05*v_['CTT']
            self.x0[ix_['X'+str(I)]] = float(v_['SCTT'])
            self.x0[ix_['T'+str(I)]] = float(v_['SCTT'])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eCOS', iet_)
        elftv = loaset(elftv,it,0,'T')
        [it,iet_,_] = s2mpj_ii( 'eSIN', iet_)
        elftv = loaset(elftv,it,0,'T')
        [it,iet_,_] = s2mpj_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'U')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['0']),int(v_['NI'])+1):
            ename = 'C'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eCOS')
            ielftype = arrset(ielftype,ie,iet_["eCOS"])
            vname = 'T'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='T')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'S'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSIN')
            ielftype = arrset(ielftype,ie,iet_["eSIN"])
            vname = 'T'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='T')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'USQ'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSQ')
            ielftype = arrset(ielftype,ie,iet_["eSQ"])
            vname = 'U'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='U')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['0']),int(v_['NI-1'])+1):
            v_['I+1'] = 1+I
            ig = ig_['EX'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['S'+str(int(v_['I+1']))])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['-H/2']))
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['S'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['-H/2']))
            ig = ig_['ENERGY']
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['USQ'+str(int(v_['I+1']))])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['H/2']))
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['USQ'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['H/2']))
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['C'+str(int(v_['I+1']))])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['AH/2']))
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['C'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['AH/2']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN(10)           345.0301196587
# LO SOLTN(50)           344.8673691861
# LO SOLTN(100)          344.8801831150
# LO SOLTN(500)          344.8748539754
# LO SOLTN(1000)         344.8788169123
# LO SOLTN(5000)         
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-COOR2-MN-V-V"
        self.objderlvl = 2
        self.conderlvl = [2]

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
        f_   = EV_[0,0]*EV_[0,0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[0,0]+EV_[0,0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eCOS(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        CC = np.cos(EV_[0,0])
        SS = np.sin(EV_[0,0])
        f_   = CC
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -SS
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = -CC
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eSIN(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        CC = np.cos(EV_[0,0])
        SS = np.sin(EV_[0,0])
        f_   = SS
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = CC
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = -SS
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

