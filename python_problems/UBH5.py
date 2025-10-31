from s2mpjlib import *
class  UBH5(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : UBH5
#    *********
# 
#    The problem is to minimize the integral of the control magnitude needed
#    to bring a vehicle, from given position and velocity, to the origin with
#    zero velocity in a fixed amount of time.  The controls are the components
#    of the vehicle acceleration. The discretization uses the trapezoidal rule.
#    This version of the problem is a variant of UBH1, where the cumulative
#    value of the objective is maintained as an additional state variable.
# 
#    The problem is convex.
# 
#    Source: unscaled problem 5 
#    (ODE = 1, CLS = 2, GRD = 1, MET = T, SEED = 0.) in
#    J.T. Betts and W.P. Huffman,
#    "Sparse Nonlinear Programming Test Problems (Release 1.0)",
#    Boeing Computer services, Seattle, July 1993.
# 
#    SIF input: Ph.L. Toint, October 1993.
# 
#    classification = "C-CLQR2-MN-V-V"
# 
#    Number of grid points
# 
#           Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER n=100, m=70    original value
# IE N                   100            $-PARAMETER n=1000, m=700
# IE N                   500            $-PARAMETER n=5000, m=3500
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'UBH5'

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
# IE N                   1000           $-PARAMETER n=10000, m=7000
# IE N                   2000           $-PARAMETER n=20000, m=14000
        v_['T0'] = 0.0
        v_['TF'] = 1000.0
        v_['RN'] = float(v_['N'])
        v_['TTIME'] = v_['TF']-v_['T0']
        v_['K'] = v_['TTIME']/v_['RN']
        v_['-K/2'] = -0.5*v_['K']
        v_['0'] = 0
        v_['1'] = 1
        v_['2'] = 2
        v_['3'] = 3
        v_['4'] = 4
        v_['5'] = 5
        v_['6'] = 6
        v_['7'] = 7
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['7'])+1):
            for T in range(int(v_['0']),int(v_['N'])+1):
                [iv,ix_,_] = s2mpj_ii('Y'+str(I)+','+str(T),ix_)
                self.xnames=arrset(self.xnames,iv,'Y'+str(I)+','+str(T))
        for I in range(int(v_['1']),int(v_['3'])+1):
            for T in range(int(v_['0']),int(v_['N'])+1):
                [iv,ix_,_] = s2mpj_ii('U'+str(I)+','+str(T),ix_)
                self.xnames=arrset(self.xnames,iv,'U'+str(I)+','+str(T))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Y'+str(int(v_['7']))+','+str(int(v_['N']))]])
        valA = np.append(valA,float(1.0))
        for I in range(int(v_['1']),int(v_['3'])+1):
            v_['I+3'] = 3+I
            for T in range(int(v_['1']),int(v_['N'])+1):
                v_['T-1'] = -1+T
                [ig,ig_,_] = s2mpj_ii('S'+str(I)+','+str(T),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'S'+str(I)+','+str(T))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['Y'+str(I)+','+str(T)]])
                valA = np.append(valA,float(1.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['Y'+str(I)+','+str(int(v_['T-1']))]])
                valA = np.append(valA,float(-1.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['Y'+str(int(v_['I+3']))+','+str(int(v_['T-1']))]])
                valA = np.append(valA,float(v_['-K/2']))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['Y'+str(int(v_['I+3']))+','+str(T)]])
                valA = np.append(valA,float(v_['-K/2']))
        for I in range(int(v_['1']),int(v_['3'])+1):
            v_['I+3'] = 3+I
            for T in range(int(v_['1']),int(v_['N'])+1):
                v_['T-1'] = -1+T
                [ig,ig_,_] = s2mpj_ii('S'+str(int(v_['I+3']))+','+str(T),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'S'+str(int(v_['I+3']))+','+str(T))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['Y'+str(int(v_['I+3']))+','+str(T)]])
                valA = np.append(valA,float(1.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['Y'+str(int(v_['I+3']))+','+str(int(v_['T-1']))]])
                valA = np.append(valA,float(-1.0))
                [ig,ig_,_] = s2mpj_ii('S'+str(int(v_['I+3']))+','+str(T),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'S'+str(int(v_['I+3']))+','+str(T))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['U'+str(I)+','+str(int(v_['T-1']))]])
                valA = np.append(valA,float(v_['-K/2']))
                [ig,ig_,_] = s2mpj_ii('S'+str(int(v_['I+3']))+','+str(T),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'S'+str(int(v_['I+3']))+','+str(T))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['U'+str(I)+','+str(T)]])
                valA = np.append(valA,float(v_['-K/2']))
        for T in range(int(v_['1']),int(v_['N'])+1):
            v_['T-1'] = -1+T
            [ig,ig_,_] = s2mpj_ii('S'+str(int(v_['7']))+','+str(T),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'S'+str(int(v_['7']))+','+str(T))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['Y'+str(int(v_['7']))+','+str(T)]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['Y'+str(int(v_['7']))+','+str(int(v_['T-1']))]])
            valA = np.append(valA,float(-1.0))
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
        for I in range(int(v_['1']),int(v_['3'])+1):
            for T in range(int(v_['0']),int(v_['N'])+1):
                self.xlower[ix_['U'+str(I)+','+str(T)]] = -1.0
                self.xupper[ix_['U'+str(I)+','+str(T)]] = 1.0
        self.xlower[ix_['Y'+str(int(v_['1']))+','+str(int(v_['0']))]] = 1000.0
        self.xupper[ix_['Y'+str(int(v_['1']))+','+str(int(v_['0']))]] = 1000.0
        self.xlower[ix_['Y'+str(int(v_['2']))+','+str(int(v_['0']))]] = 1000.0
        self.xupper[ix_['Y'+str(int(v_['2']))+','+str(int(v_['0']))]] = 1000.0
        self.xlower[ix_['Y'+str(int(v_['3']))+','+str(int(v_['0']))]] = 1000.0
        self.xupper[ix_['Y'+str(int(v_['3']))+','+str(int(v_['0']))]] = 1000.0
        self.xlower[ix_['Y'+str(int(v_['4']))+','+str(int(v_['0']))]] = -10.0
        self.xupper[ix_['Y'+str(int(v_['4']))+','+str(int(v_['0']))]] = -10.0
        self.xlower[ix_['Y'+str(int(v_['5']))+','+str(int(v_['0']))]] = 10.0
        self.xupper[ix_['Y'+str(int(v_['5']))+','+str(int(v_['0']))]] = 10.0
        self.xlower[ix_['Y'+str(int(v_['6']))+','+str(int(v_['0']))]] = -10.0
        self.xupper[ix_['Y'+str(int(v_['6']))+','+str(int(v_['0']))]] = -10.0
        self.xlower[ix_['Y'+str(int(v_['7']))+','+str(int(v_['0']))]] = 0.0
        self.xupper[ix_['Y'+str(int(v_['7']))+','+str(int(v_['0']))]] = 0.0
        for I in range(int(v_['1']),int(v_['6'])+1):
            self.xlower[ix_['Y'+str(I)+','+str(int(v_['N']))]] = 0.0
            self.xupper[ix_['Y'+str(I)+','+str(int(v_['N']))]] = 0.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'V')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for T in range(int(v_['0']),int(v_['N'])+1):
            for I in range(int(v_['1']),int(v_['3'])+1):
                ename = 'E'+str(I)+','+str(T)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eSQ')
                ielftype = arrset(ielftype,ie,iet_["eSQ"])
                vname = 'U'+str(I)+','+str(T)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='V')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for T in range(int(v_['1']),int(v_['N'])+1):
            v_['T-1'] = -1+T
            for I in range(int(v_['1']),int(v_['3'])+1):
                ig = ig_['S'+str(int(v_['7']))+','+str(T)]
                posel = len(self.grelt[ig])
                self.grelt  = (
                      loaset(self.grelt,ig,posel,ie_['E'+str(I)+','+str(int(v_['T-1']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(v_['-K/2']))
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(I)+','+str(T)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(v_['-K/2']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 0.0
#    Solution
# LO SOLTN(10)           1.14735202967
# LO SOLTN(100)          1.11631518169
# LO SOLTN(1000)         1.11598643493
# LO SOLTN(2000)         1.11587382445
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
        self.pbclass   = "C-CLQR2-MN-V-V"
        self.x0        = np.zeros((self.n,1))
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
        f_   = EV_[0]*EV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[0]+EV_[0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

