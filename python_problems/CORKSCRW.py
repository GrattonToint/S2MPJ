from s2mpjlib import *
class  CORKSCRW(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CORKSCRW
#    *********
# 
#    A nonlinear optimal control problem with both state- and
#    control constraints.
#    The problem is to control (using an applied force of limited
#    magnitude) a mass moving in the 3D space, such that its
#    trajectory lies within a prescribed distance TOL of the
#    corkscreww-like curve defined by
#               y = sin(x), z = cos(x),
#    and such that it stops at a given abscissa XT in minimum time.
#    The mass is initially stationary at (0,0,1).
# 
#    Source:
#    Ph. Toint, private communication.
# 
#    SIF input: Ph. Toint, April 1991.
# 
#    classification = "C-CSOR2-AN-V-V"
# 
#    Number of time intervals
#    The number of variables is 9T+6, of which 9 are fixed.
# 
#           Alternative values for the SIF file parameters:
# IE T                   10             $-PARAMETER n = 96     original value
# IE T                   50             $-PARAMETER n = 456
# IE T                   100            $-PARAMETER n = 906
# IE T                   500            $-PARAMETER n = 4506
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'CORKSCRW'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['T'] = int(10);  #  SIF file default value
        else:
            v_['T'] = int(args[0])
# IE T                   1000           $-PARAMETER n = 9006
        if nargin<2:
            v_['XT'] = float(10.0);  #  SIF file default value
        else:
            v_['XT'] = float(args[1])
        if nargin<3:
            v_['MASS'] = float(0.37);  #  SIF file default value
        else:
            v_['MASS'] = float(args[2])
        if nargin<4:
            v_['TOL'] = float(0.1);  #  SIF file default value
        else:
            v_['TOL'] = float(args[3])
        v_['0'] = 0
        v_['1'] = 1
        v_['RT'] = float(v_['T'])
        v_['T+1'] = 1.0+v_['RT']
        v_['H'] = v_['XT']/v_['RT']
        v_['1/H'] = 1.0/v_['H']
        v_['-1/H'] = -1.0*v_['1/H']
        v_['M/H'] = v_['MASS']/v_['H']
        v_['-M/H'] = -1.0*v_['M/H']
        v_['TOLSQ'] = v_['TOL']*v_['TOL']
        v_['XTT+1'] = v_['XT']*v_['T+1']
        v_['W'] = 0.5*v_['XTT+1']
        for I in range(int(v_['1']),int(v_['T'])+1):
            v_['RI'] = float(I)
            v_['TI'] = v_['RI']*v_['H']
            v_['W/T'+str(I)] = v_['W']/v_['TI']
        v_['FMAX'] = v_['XT']/v_['RT']
        v_['-FMAX'] = -1.0*v_['FMAX']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['0']),int(v_['T'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I))
            [iv,ix_,_] = s2mpj_ii('Y'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'Y'+str(I))
            [iv,ix_,_] = s2mpj_ii('Z'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'Z'+str(I))
            [iv,ix_,_] = s2mpj_ii('VX'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'VX'+str(I))
            [iv,ix_,_] = s2mpj_ii('VY'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'VY'+str(I))
            [iv,ix_,_] = s2mpj_ii('VZ'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'VZ'+str(I))
        for I in range(int(v_['1']),int(v_['T'])+1):
            [iv,ix_,_] = s2mpj_ii('UX'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'UX'+str(I))
            [iv,ix_,_] = s2mpj_ii('UY'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'UY'+str(I))
            [iv,ix_,_] = s2mpj_ii('UZ'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'UZ'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['1']),int(v_['T'])+1):
            [ig,ig_,_] = s2mpj_ii('OX'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            self.gscale = arrset(self.gscale,ig,float(v_['W/T'+str(I)]))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(I)]])
            valA = np.append(valA,float(1.0))
        for I in range(int(v_['1']),int(v_['T'])+1):
            v_['I-1'] = -1+I
            [ig,ig_,_] = s2mpj_ii('ACX'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'ACX'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['VX'+str(I)]])
            valA = np.append(valA,float(v_['M/H']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['VX'+str(int(v_['I-1']))]])
            valA = np.append(valA,float(v_['-M/H']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['UX'+str(I)]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('ACY'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'ACY'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['VY'+str(I)]])
            valA = np.append(valA,float(v_['M/H']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['VY'+str(int(v_['I-1']))]])
            valA = np.append(valA,float(v_['-M/H']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['UY'+str(I)]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('ACZ'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'ACZ'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['VZ'+str(I)]])
            valA = np.append(valA,float(v_['M/H']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['VZ'+str(int(v_['I-1']))]])
            valA = np.append(valA,float(v_['-M/H']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['UZ'+str(I)]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('PSX'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'PSX'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(I)]])
            valA = np.append(valA,float(v_['1/H']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['I-1']))]])
            valA = np.append(valA,float(v_['-1/H']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['VX'+str(I)]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('PSY'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'PSY'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['Y'+str(I)]])
            valA = np.append(valA,float(v_['1/H']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['Y'+str(int(v_['I-1']))]])
            valA = np.append(valA,float(v_['-1/H']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['VY'+str(I)]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('PSZ'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'PSZ'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['Z'+str(I)]])
            valA = np.append(valA,float(v_['1/H']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['Z'+str(int(v_['I-1']))]])
            valA = np.append(valA,float(v_['-1/H']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['VZ'+str(I)]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('SC'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'SC'+str(I))
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
        for I in range(int(v_['1']),int(v_['T'])+1):
            self.gconst = arrset(self.gconst,ig_['OX'+str(I)],float(v_['XT']))
            self.gconst = arrset(self.gconst,ig_['SC'+str(I)],float(v_['TOLSQ']))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        self.xlower[ix_['X'+str(int(v_['0']))]] = 0.0
        self.xupper[ix_['X'+str(int(v_['0']))]] = 0.0
        self.xlower[ix_['Y'+str(int(v_['0']))]] = 0.0
        self.xupper[ix_['Y'+str(int(v_['0']))]] = 0.0
        self.xlower[ix_['Z'+str(int(v_['0']))]] = 1.0
        self.xupper[ix_['Z'+str(int(v_['0']))]] = 1.0
        self.xlower[ix_['VX'+str(int(v_['0']))]] = 0.0
        self.xupper[ix_['VX'+str(int(v_['0']))]] = 0.0
        self.xlower[ix_['VY'+str(int(v_['0']))]] = 0.0
        self.xupper[ix_['VY'+str(int(v_['0']))]] = 0.0
        self.xlower[ix_['VZ'+str(int(v_['0']))]] = 0.0
        self.xupper[ix_['VZ'+str(int(v_['0']))]] = 0.0
        self.xlower[ix_['VX'+str(int(v_['T']))]] = 0.0
        self.xupper[ix_['VX'+str(int(v_['T']))]] = 0.0
        self.xlower[ix_['VY'+str(int(v_['T']))]] = 0.0
        self.xupper[ix_['VY'+str(int(v_['T']))]] = 0.0
        self.xlower[ix_['VZ'+str(int(v_['T']))]] = 0.0
        self.xupper[ix_['VZ'+str(int(v_['T']))]] = 0.0
        for I in range(int(v_['1']),int(v_['T'])+1):
            self.xlower[ix_['UX'+str(I)]] = v_['-FMAX']
            self.xupper[ix_['UX'+str(I)]] = v_['FMAX']
            self.xlower[ix_['UY'+str(I)]] = v_['-FMAX']
            self.xupper[ix_['UY'+str(I)]] = v_['FMAX']
            self.xlower[ix_['UZ'+str(I)]] = v_['-FMAX']
            self.xupper[ix_['UZ'+str(I)]] = v_['FMAX']
            self.xlower[ix_['X'+str(I)]] = 0.0
            self.xupper[ix_['X'+str(I)]] = v_['XT']
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        self.x0[ix_['X'+str(int(v_['0']))]] = float(0.0)
        self.x0[ix_['Y'+str(int(v_['0']))]] = float(0.0)
        self.x0[ix_['Z'+str(int(v_['0']))]] = float(1.0)
        self.x0[ix_['VX'+str(int(v_['0']))]] = float(0.0)
        self.x0[ix_['VY'+str(int(v_['0']))]] = float(0.0)
        self.x0[ix_['VZ'+str(int(v_['0']))]] = float(0.0)
        self.x0[ix_['VX'+str(int(v_['T']))]] = float(0.0)
        self.x0[ix_['VY'+str(int(v_['T']))]] = float(0.0)
        self.x0[ix_['VZ'+str(int(v_['T']))]] = float(0.0)
        for I in range(int(v_['1']),int(v_['T'])+1):
            v_['RI'] = float(I)
            v_['TI'] = v_['RI']*v_['H']
            if('X'+str(I) in ix_):
                self.x0[ix_['X'+str(I)]] = float(v_['TI'])
            else:
                self.y0  = (
                      arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X'+str(I)]),float(v_['TI'])))
            if('VX'+str(I) in ix_):
                self.x0[ix_['VX'+str(I)]] = float(1.0)
            else:
                self.y0  = (
                      arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['VX'+str(I)]),float(1.0)))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eERRSIN', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        [it,iet_,_] = s2mpj_ii( 'eERRCOS', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Z')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['1']),int(v_['T'])+1):
            ename = 'ES'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eERRSIN')
            ielftype = arrset(ielftype,ie,iet_["eERRSIN"])
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'Y'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='Y')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'EC'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eERRCOS')
            ielftype = arrset(ielftype,ie,iet_["eERRCOS"])
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'Z'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='Z')[0]
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
        for I in range(int(v_['1']),int(v_['T'])+1):
            ig = ig_['OX'+str(I)]
            self.grftype = arrset(self.grftype,ig,'gL2')
            ig = ig_['SC'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['ES'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = posel+1
            self.grelt = loaset(self.grelt,ig,posel,ie_['EC'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel, 1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN(10)           1.1601050195
# LO SOLTN(50)           26.484181830
# LO SOLTN(100)          44.368110588
# LO SOLTN(500)
# LO SOLTN(1000)
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.cupper[np.arange(self.nle)] = np.zeros((self.nle,1))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CSOR2-AN-V-V"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eERRSIN(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        SINX = np.sin(EV_[0])
        COSX = np.cos(EV_[0])
        ERR = EV_[1]-SINX
        f_   = ERR*ERR
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -2.0*ERR*COSX
            g_[1] = 2.0*ERR
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 2.0*(COSX**2+ERR*SINX)
                H_[0,1] = -2.0*COSX
                H_[1,0] = H_[0,1]
                H_[1,1] = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eERRCOS(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        SINX = np.sin(EV_[0])
        COSX = np.cos(EV_[0])
        ERR = EV_[1]-COSX
        f_   = ERR*ERR
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*ERR*SINX
            g_[1] = 2.0*ERR
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 2.0*(SINX**2+ERR*COSX)
                H_[0,1] = 2.0*SINX
                H_[1,0] = H_[0,1]
                H_[1,1] = 2.0
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

