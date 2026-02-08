from s2mpjlib import *
class  DRUGDISE(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DRUGDISE
#    *********
# 
#    This is a variant of the drug displacement problem DRUGDIS where the
#    state equations have been Expanded in term of more intermediate
#    functions, each one of them being less nonlinear.
# 
#    The problem is based on the kinetic model of Aarons and Rowland which
#    simulates the interaction of the two drugs (warfarin and phenylnutazone)
#    in a patient bloodstream.  The state variable are the concentrations of
#    unbound warfarin (w) and phenylbutazone (p).  The problem is to control
#    the rate of injection (u) of the pain-killing phenylbutazone so that both
#    drugs reach a specified steady-state in minimum time and the concentration
#    of warfarin does not rise above a toxicity level.
# 
#    The problem is discretized using the trapeziodal rule.  It is non-convex.
# 
#    Source:
#    H. Maurer and M. Wiegand,
#    "Numerical solution of a drug displacement problem with bounded state
#    variables",
#    Optimal Control Applications and Methods 13, pp. 43-55, 1992.
# 
#    SIF input: Ph. Toint, Nov 1993.
# 
#    classification = "C-CLOR2-MY-V-V"
# 
#    Discretization: specify the number of interior points + 1
# 
#           Alternative values for the SIF file parameters:
# IE NI                  10             $-PARAMETER n=63, m=50 
# IE NI                  100            $-PARAMETER n=603, m=500   original value
# IE NI                  100            $-PARAMETER n=6003, m=5000 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'DRUGDISE'

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
        if nargin<2:
            v_['TOXIC'] = float(0.026);  #  SIF file default value
        else:
            v_['TOXIC'] = float(args[1])
        if nargin<3:
            v_['WSS'] = float(0.02);  #  SIF file default value
        else:
            v_['WSS'] = float(args[2])
        if nargin<4:
            v_['UMAX'] = float(8.0);  #  SIF file default value
        else:
            v_['UMAX'] = float(args[3])
        if nargin<5:
            v_['PSTART'] = float(0.0);  #  SIF file default value
        else:
            v_['PSTART'] = float(args[4])
        if nargin<6:
            v_['PFINAL'] = float(2.0);  #  SIF file default value
        else:
            v_['PFINAL'] = float(args[5])
        if nargin<7:
            v_['Z'] = float(46.4);  #  SIF file default value
        else:
            v_['Z'] = float(args[6])
        v_['AVP'] = v_['PSTART']+v_['PFINAL']
        v_['AVP'] = 0.5*v_['AVP']
        v_['-Z'] = -1.0*v_['Z']
        v_['-ZZ'] = v_['Z']*v_['-Z']
        v_['NI-1'] = -1+v_['NI']
        v_['RNI'] = float(v_['NI'])
        v_['-1/NI'] = -1.0/v_['RNI']
        v_['-Z/NI'] = v_['Z']*v_['-1/NI']
        v_['0'] = 0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [iv,ix_,_] = s2mpj_ii('TF',ix_)
        self.xnames=arrset(self.xnames,iv,'TF')
        self.xscale = arrset(self.xscale,iv,200.0)
        for I in range(int(v_['0']),int(v_['NI'])+1):
            [iv,ix_,_] = s2mpj_ii('W'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'W'+str(I))
            self.xscale = arrset(self.xscale,iv,0.02)
        for I in range(int(v_['0']),int(v_['NI'])+1):
            [iv,ix_,_] = s2mpj_ii('P'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'P'+str(I))
        for I in range(int(v_['0']),int(v_['NI-1'])+1):
            [iv,ix_,_] = s2mpj_ii('U'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'U'+str(I))
        for I in range(int(v_['0']),int(v_['NI-1'])+1):
            [iv,ix_,_] = s2mpj_ii('A'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'A'+str(I))
            self.xscale = arrset(self.xscale,iv,200.0)
        for I in range(int(v_['0']),int(v_['NI-1'])+1):
            [iv,ix_,_] = s2mpj_ii('B'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'B'+str(I))
            self.xscale = arrset(self.xscale,iv,200.0)
        for I in range(int(v_['0']),int(v_['NI-1'])+1):
            [iv,ix_,_] = s2mpj_ii('C'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'C'+str(I))
            self.xscale = arrset(self.xscale,iv,0.0000001)
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('TFINAL',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['TF']])
        valA = np.append(valA,float(1.0))
        self.gscale = arrset(self.gscale,ig,float(100.0))
        for I in range(int(v_['0']),int(v_['NI-1'])+1):
            v_['I+1'] = 1+I
            [ig,ig_,_] = s2mpj_ii('EW'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'EW'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['W'+str(int(v_['I+1']))]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['W'+str(I)]])
            valA = np.append(valA,float(-1.0))
            self.gscale = arrset(self.gscale,ig,float(0.02))
            [ig,ig_,_] = s2mpj_ii('EP'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'EP'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['P'+str(int(v_['I+1']))]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['P'+str(I)]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('EA'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'EA'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['A'+str(I)]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['P'+str(I)]])
            valA = np.append(valA,float(v_['-Z']))
            self.gscale = arrset(self.gscale,ig,float(200.0))
            [ig,ig_,_] = s2mpj_ii('EB'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'EB'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['B'+str(I)]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['W'+str(I)]])
            valA = np.append(valA,float(v_['-Z']))
            self.gscale = arrset(self.gscale,ig,float(200.0))
            [ig,ig_,_] = s2mpj_ii('EC'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'EC'+str(I))
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
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        for I in range(int(v_['0']),int(v_['NI-1'])+1):
            self.gconst = arrset(self.gconst,ig_['EA'+str(I)],float(232.0))
            self.gconst = arrset(self.gconst,ig_['EB'+str(I)],float(232.0))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        for I in range(int(v_['0']),int(v_['NI-1'])+1):
            self.xlower[ix_['C'+str(I)]] = -float('Inf')
            self.xupper[ix_['C'+str(I)]] = +float('Inf')
        self.xlower[ix_['TF']] = 200.0
        for I in range(int(v_['0']),int(v_['NI'])+1):
            self.xupper[ix_['W'+str(I)]] = v_['TOXIC']
        for I in range(int(v_['0']),int(v_['NI-1'])+1):
            self.xupper[ix_['U'+str(I)]] = v_['UMAX']
        self.xlower[ix_['W'+str(int(v_['0']))]] = v_['WSS']
        self.xupper[ix_['W'+str(int(v_['0']))]] = v_['WSS']
        self.xlower[ix_['W'+str(int(v_['NI']))]] = v_['WSS']
        self.xupper[ix_['W'+str(int(v_['NI']))]] = v_['WSS']
        self.xlower[ix_['P'+str(int(v_['0']))]] = v_['PSTART']
        self.xupper[ix_['P'+str(int(v_['0']))]] = v_['PSTART']
        self.xlower[ix_['P'+str(int(v_['NI']))]] = v_['PFINAL']
        self.xupper[ix_['P'+str(int(v_['NI']))]] = v_['PFINAL']
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        v_['2W/10'] = 0.2*v_['WSS']
        v_['2P/10'] = 0.2*v_['AVP']
        v_['2(W+P)/10'] = v_['2W/10']+v_['2P/10']
        v_['D'] = 1.0+v_['2(W+P)/10']
        v_['DD'] = v_['D']*v_['D']
        v_['ZP'] = v_['AVP']*v_['Z']
        v_['ZW'] = v_['WSS']*v_['Z']
        v_['AA'] = v_['DD']+v_['ZP']
        v_['AA'] = 232.0+v_['AA']
        v_['BB'] = v_['DD']+v_['ZW']
        v_['BB'] = 232.0+v_['BB']
        v_['AB'] = v_['AA']*v_['BB']
        v_['WP'] = v_['WSS']*v_['AVP']
        v_['-ZZWP'] = v_['WP']*v_['-ZZ']
        v_['CD'] = v_['AB']+v_['-ZZWP']
        v_['CC'] = v_['DD']/v_['CD']
        for I in range(int(v_['0']),int(v_['NI-1'])+1):
            self.x0[ix_['W'+str(I)]] = float(v_['WSS'])
            self.x0[ix_['P'+str(I)]] = float(v_['AVP'])
            self.x0[ix_['U'+str(I)]] = float(v_['UMAX'])
            self.x0[ix_['A'+str(I)]] = float(v_['AA'])
            self.x0[ix_['B'+str(I)]] = float(v_['BB'])
            self.x0[ix_['C'+str(I)]] = float(v_['CC'])
        self.x0[ix_['TF']] = float(240.0)
        self.x0[ix_['W'+str(int(v_['NI']))]] = float(v_['WSS'])
        self.x0[ix_['P'+str(int(v_['NI']))]] = float(v_['PFINAL'])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'en3S', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftv = loaset(elftv,it,3,'V4')
        [it,iet_,_] = s2mpj_ii( 'en3D2', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftv = loaset(elftv,it,3,'V4')
        elftv = loaset(elftv,it,4,'V5')
        [it,iet_,_] = s2mpj_ii( 'eDSQ', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        [it,iet_,_] = s2mpj_ii( 'en3PR', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['0']),int(v_['NI-1'])+1):
            ename = 'WA'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'en3S')
            ielftype = arrset(ielftype,ie,iet_["en3S"])
            vname = 'TF'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'C'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'A'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'W'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V4')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'WB'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'en3D2')
            ielftype = arrset(ielftype,ie,iet_["en3D2"])
            vname = 'TF'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'C'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'W'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'U'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V4')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'P'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V5')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'PA'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'en3D2')
            ielftype = arrset(ielftype,ie,iet_["en3D2"])
            vname = 'TF'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'C'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'B'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'U'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V4')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'P'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V5')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'PB'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'en3S')
            ielftype = arrset(ielftype,ie,iet_["en3S"])
            vname = 'TF'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'C'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'P'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'W'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V4')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'DD'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eDSQ')
            ielftype = arrset(ielftype,ie,iet_["eDSQ"])
            vname = 'W'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'P'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'CA'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'en3PR')
            ielftype = arrset(ielftype,ie,iet_["en3PR"])
            vname = 'C'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'A'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'B'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'CB'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'en3PR')
            ielftype = arrset(ielftype,ie,iet_["en3PR"])
            vname = 'C'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'P'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'W'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['0']),int(v_['NI-1'])+1):
            ig = ig_['EW'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['WA'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['-1/NI']))
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['WB'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['-Z/NI']))
            ig = ig_['EP'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['PA'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['-1/NI']))
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['PB'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['-Z/NI']))
            ig = ig_['EA'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['DD'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
            ig = ig_['EB'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['DD'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
            ig = ig_['EC'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['CA'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(1.0))
            posel = posel+1
            self.grelt = loaset(self.grelt,ig,posel,ie_['DD'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['CB'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['-ZZ']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 200.0
#    Solution
# LO SOLTN               ????
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
        self.pbclass   = "C-CLOR2-MY-V-V"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def e_globs(self):

        import numpy as np
        self.efpar = np.array([]);
        self.efpar = arrset( self.efpar,0,0.02)
        return pbm

    @staticmethod
    def en3S(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        WSSMV4 = self.efpar[0]-EV_[3,0]
        f_   = EV_[0,0]*EV_[1,0]*EV_[2,0]*WSSMV4
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1,0]*EV_[2,0]*WSSMV4
            g_[1] = EV_[0,0]*EV_[2,0]*WSSMV4
            g_[2] = EV_[0,0]*EV_[1,0]*WSSMV4
            g_[3] = -EV_[0,0]*EV_[1,0]*EV_[2,0]
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,1] = EV_[2,0]*WSSMV4
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1,0]*WSSMV4
                H_[2,0] = H_[0,2]
                H_[0,3] = -EV_[1,0]*EV_[2,0]
                H_[3,0] = H_[0,3]
                H_[1,2] = EV_[0,0]*WSSMV4
                H_[2,1] = H_[1,2]
                H_[1,3] = -EV_[0,0]*EV_[2,0]
                H_[3,1] = H_[1,3]
                H_[2,3] = -EV_[0,0]*EV_[1,0]
                H_[3,2] = H_[2,3]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def en3D2(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((4,5))
        IV_ = np.zeros(4)
        U_[0,0] = U_[0,0]+1
        U_[1,1] = U_[1,1]+1
        U_[2,2] = U_[2,2]+1
        U_[3,3] = U_[3,3]+1
        U_[3,4] = U_[3,4]-2
        IV_[0] = to_scalar(U_[0:1,:].dot(EV_))
        IV_[1] = to_scalar(U_[1:2,:].dot(EV_))
        IV_[2] = to_scalar(U_[2:3,:].dot(EV_))
        IV_[3] = to_scalar(U_[3:4,:].dot(EV_))
        f_   = IV_[0]*IV_[1]*IV_[2]*IV_[3]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = IV_[1]*IV_[2]*IV_[3]
            g_[1] = IV_[0]*IV_[2]*IV_[3]
            g_[2] = IV_[0]*IV_[1]*IV_[3]
            g_[3] = IV_[0]*IV_[1]*IV_[2]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,1] = IV_[2]*IV_[3]
                H_[1,0] = H_[0,1]
                H_[0,2] = IV_[1]*IV_[3]
                H_[2,0] = H_[0,2]
                H_[0,3] = IV_[1]*IV_[2]
                H_[3,0] = H_[0,3]
                H_[1,2] = IV_[0]*IV_[3]
                H_[2,1] = H_[1,2]
                H_[1,3] = IV_[0]*IV_[2]
                H_[3,1] = H_[1,3]
                H_[2,3] = IV_[0]*IV_[1]
                H_[3,2] = H_[2,3]
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eDSQ(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((1,2))
        IV_ = np.zeros(1)
        U_[0,0] = U_[0,0]+2.000000e-01
        U_[0,1] = U_[0,1]+2.000000e-01
        IV_[0] = to_scalar(U_[0:1,:].dot(EV_))
        f_   = IV_[0]*IV_[0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = IV_[0]+IV_[0]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def en3PR(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0,0]*EV_[1,0]*EV_[2,0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1,0]*EV_[2,0]
            g_[1] = EV_[0,0]*EV_[2,0]
            g_[2] = EV_[0,0]*EV_[1,0]
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = EV_[2,0]
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1,0]
                H_[2,0] = H_[0,2]
                H_[1,2] = EV_[0,0]
                H_[2,1] = H_[1,2]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

