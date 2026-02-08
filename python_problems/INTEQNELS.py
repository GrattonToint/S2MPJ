from s2mpjlib import *
class  INTEQNELS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : INTEQNELS
#    *********
#    The discrete integral problem (INTEGREQ) without fixed variables
#    in least-squares form.
# 
#    Source:  Problem 29 in
#    J.J. More', B.S. Garbow and K.E. Hillstrom,
#    "Testing Unconstrained Optimization Software",
#    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
# 
#    SIF input: Ph. Toint, Feb 1990.
#    Modification to remove fixed variables: Nick Gould, Oct 2015.
# 
#    classification = "C-CSUR2-AN-V-0"
# 
#    N+2 is the number of discretization points .
#    The number of free variables is N.
# 
#           Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER
# IE N                   50             $-PARAMETER     original value
# IE N                   100            $-PARAMETER
# IE N                   500            $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'INTEQNELS'

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
        v_['0'] = 0
        v_['1'] = 1
        v_['N+1'] = 1+v_['N']
        v_['RN+1'] = float(v_['N+1'])
        v_['H'] = 1.0/v_['RN+1']
        v_['HALFH'] = 0.5e0*v_['H']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['0']),int(v_['N+1'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('G'+str(int(v_['0'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X'+str(int(v_['0']))]])
        valA = np.append(valA,float(1.0))
        for I in range(int(v_['1']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('G'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(I)]])
            valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('G'+str(int(v_['N+1'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X'+str(int(v_['N+1']))]])
        valA = np.append(valA,float(1.0))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        selfnob      = ngrp
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.x0[ix_['X'+str(int(v_['0']))]] = float(0.0)
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['REALI'] = float(I)
            v_['IH'] = v_['REALI']*v_['H']
            v_['IH-1'] = -1.0+v_['IH']
            v_['TI'] = v_['IH']*v_['IH-1']
            self.x0[ix_['X'+str(I)]] = float(v_['TI'])
        self.x0[ix_['X'+str(int(v_['N+1']))]] = float(0.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eVBCUBE', iet_)
        elftv = loaset(elftv,it,0,'V')
        elftp = []
        elftp = loaset(elftp,it,0,'B')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for J in range(int(v_['1']),int(v_['N'])+1):
            v_['REALJ'] = float(J)
            v_['TJ'] = v_['REALJ']*v_['H']
            v_['1+TJ'] = 1.0+v_['TJ']
            ename = 'A'+str(J)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eVBCUBE')
            ielftype = arrset(ielftype,ie,iet_["eVBCUBE"])
            vname = 'X'+str(J)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='B')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['1+TJ']))
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
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['REALI'] = float(I)
            v_['TI'] = v_['REALI']*v_['H']
            v_['-TI'] = -1.0*v_['TI']
            v_['1-TI'] = 1.0+v_['-TI']
            v_['P1'] = v_['1-TI']*v_['HALFH']
            for J in range(int(v_['1']),int(I)+1):
                v_['REALJ'] = float(J)
                v_['TJ'] = v_['REALJ']*v_['H']
                v_['WIL'] = v_['P1']*v_['TJ']
                ig = ig_['G'+str(I)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['A'+str(J)])
                self.grelw = loaset(self.grelw,ig,posel,float(v_['WIL']))
            v_['I+1'] = 1+I
            v_['P2'] = v_['TI']*v_['HALFH']
            for J in range(int(v_['I+1']),int(v_['N'])+1):
                v_['REALJ'] = float(J)
                v_['TJ'] = v_['REALJ']*v_['H']
                v_['-TJ'] = -1.0*v_['TJ']
                v_['1-TJ'] = 1.0+v_['-TJ']
                v_['WIU'] = v_['P2']*v_['1-TJ']
                ig = ig_['G'+str(I)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['A'+str(J)])
                self.grelw = loaset(self.grelw,ig,posel,float(v_['WIU']))
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass   = "C-CSUR2-AN-V-0"
        self.objderlvl = 2


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eVBCUBE(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        VPLUSB = EV_[0,0]+self.elpar[iel_][0]
        f_   = VPLUSB**3
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 3.0*VPLUSB**2
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 6.0*VPLUSB
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

