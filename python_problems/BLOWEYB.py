from s2mpjlib import *
class  BLOWEYB(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : BLOWEYB
#    *********
# 
#    A nonconvex quadratic program proposed by 
#    James Blowey (University of Durham)
# 
#    Given function v(s) and f(s) = v(s) + A(inv) v(s), s in [0,1],
#    minimize 
# 
#         (u(s) - v(s))(trans) ( A + A(inv) ) u(s) - (u(s) - v(s))(trans)f(s)
# 
#    where 
# 
#       u(s) in [-1,1] and int[0,1] u(s) ds = int[0,1] v(s) ds
# 
#    and A is the 
# 
#       "- Laplacian with Neumann boundary conditions on a uniform mesh"
# 
#    The troublesome term A(inv) u(s) is replaced by the additional 
#    variable w(s) and the constraint A w(s) = u(s)
# 
#    The function v(s) is chosen to be 
# 
#      0  a  b   c  d   1
#   1  ----         ----
#          \       /
#           \     /
#            \   /
#  -1         --- 
# 
#    Thus the problem is formulated as the nonconvex QP
# 
#    minimize 
# 
#         u(s) (trans) A u(s) + u(s) (trans) w(s) - 
#         v(s)(trans) A u(s) - 2.0 v(s)(trans) w(s) - 
#         u(s)(trans) v(s) + constant (ignored)
# 
#    subject to A w(s) = u(s), 
#               u(s) in [-1,1],
#           and int[0,1] u(s) ds = 1 + a + b - c - d
# 
#    Case B: a = 0.1, b = 0.4, c = 0.6 and d = 0.9. 
# 
#    Source: a simplification of
#    J.F. Blowey and C.M. Elliott,
#    "The Cahn-Hilliard gradient theory for phase separation with 
#    non-smooth free energy Part II: Numerical analysis",
#    European Journal of Applied Mathematics (3) pp 147-179, 1992.
# 
#    SIF input: Nick Gould, August 1996
# 
#    classification = "C-CQLR2-MN-V-V"
# 
#    The number of discretization intervals
# 
#           Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER  n = 22, m = 12
# IE N                   100            $-PARAMETER  n = 202, m = 102
# IE N                   1000           $-PARAMETER  n = 2002, m = 1002
# IE N                   2000           $-PARAMETER  n = 4002, m = 2002
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'BLOWEYB'

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
# IE N                   4000           $-PARAMETER  n = 8002, m = 4002
# IE N                   8000           $-PARAMETER  n = 16002, m = 8002
        v_['0'] = 0
        v_['1'] = 1
        v_['2'] = 2
        v_['3'] = 3
        v_['5'] = 5
        v_['9'] = 9
        v_['10'] = 10
        v_['ONE'] = 1.0
        v_['-ONE'] = -1.0
        v_['TWO'] = 2.0
        v_['-TWO'] = -2.0
        v_['RN'] = float(v_['N'])
        v_['N**2'] = v_['RN']*v_['RN']
        v_['N-1'] = -1+v_['N']
        v_['1/N**2'] = v_['ONE']/v_['N**2']
        v_['-1/N**2'] = v_['-ONE']/v_['N**2']
        v_['-2/N**2'] = 2.0*v_['-1/N**2']
        v_['2N**2'] = 2.0*v_['N**2']
        v_['-2N**2'] = -2.0*v_['N**2']
        v_['N/10'] = int(np.fix(v_['N']/v_['10']))
        v_['N/5'] = int(np.fix(v_['N']/v_['5']))
        v_['NA'] = v_['N/10']
        v_['A'] = float(v_['NA'])
        v_['A'] = v_['A']/v_['RN']
        v_['NA+1'] = 1+v_['NA']
        v_['NB'] = v_['N/5']*v_['2']
        v_['B'] = float(v_['NB'])
        v_['B'] = v_['B']/v_['RN']
        v_['NB+1'] = 1+v_['NB']
        v_['NC'] = v_['N/5']*v_['3']
        v_['C'] = float(v_['NC'])
        v_['C'] = v_['C']/v_['RN']
        v_['NC+1'] = 1+v_['NC']
        v_['ND'] = v_['N/10']*v_['9']
        v_['D'] = float(v_['ND'])
        v_['D'] = v_['D']/v_['RN']
        v_['ND+1'] = 1+v_['ND']
        v_['INT'] = v_['ONE']
        v_['INT'] = v_['INT']+v_['A']
        v_['INT'] = v_['INT']+v_['B']
        v_['INT'] = v_['INT']-v_['C']
        v_['INT'] = v_['INT']-v_['D']
        v_['INT'] = v_['INT']*v_['RN']
        for I in range(int(v_['0']),int(v_['NA'])+1):
            v_['V'+str(I)] = 1.0
        v_['STEP'] = v_['B']-v_['A']
        v_['STEP'] = v_['STEP']*v_['RN']
        v_['STEP'] = v_['TWO']/v_['STEP']
        for I in range(int(v_['NA+1']),int(v_['NB'])+1):
            v_['J'] = I-v_['NA']
            v_['RJ'] = float(v_['J'])
            v_['VAL'] = v_['RJ']*v_['STEP']
            v_['VAL'] = v_['ONE']-v_['VAL']
            v_['V'+str(I)] = v_['VAL']
        for I in range(int(v_['NB+1']),int(v_['NC'])+1):
            v_['V'+str(I)] = -1.0
        v_['STEP'] = v_['D']-v_['C']
        v_['STEP'] = v_['STEP']*v_['RN']
        v_['STEP'] = v_['TWO']/v_['STEP']
        for I in range(int(v_['NC+1']),int(v_['ND'])+1):
            v_['J'] = I-v_['NC']
            v_['RJ'] = float(v_['J'])
            v_['VAL'] = v_['RJ']*v_['STEP']
            v_['VAL'] = v_['-ONE']+v_['VAL']
            v_['V'+str(I)] = v_['VAL']
        for I in range(int(v_['ND']),int(v_['N'])+1):
            v_['V'+str(I)] = 1.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['0']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('U'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'U'+str(I))
            [iv,ix_,_] = s2mpj_ii('W'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'W'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['0']),int(v_['N'])+1):
            v_['VAL'] = v_['V'+str(I)]*v_['-1/N**2']
            [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['U'+str(I)]])
            valA = np.append(valA,float(v_['VAL']))
            v_['VAL'] = v_['V'+str(I)]*v_['-2/N**2']
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['W'+str(I)]])
            valA = np.append(valA,float(v_['VAL']))
        v_['VAL'] = v_['V'+str(int(v_['1']))]-v_['V'+str(int(v_['0']))]
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['U'+str(int(v_['0']))]])
        valA = np.append(valA,float(v_['VAL']))
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            v_['I+1'] = 1+I
            v_['I-1'] = -1+I
            v_['VAL'] = -2.0*v_['V'+str(I)]
            v_['VAL'] = v_['VAL']+v_['V'+str(int(v_['I-1']))]
            v_['VAL'] = v_['VAL']+v_['V'+str(int(v_['I+1']))]
            [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['U'+str(I)]])
            valA = np.append(valA,float(v_['VAL']))
        v_['VAL'] = v_['V'+str(int(v_['N-1']))]-v_['V'+str(int(v_['N']))]
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['U'+str(int(v_['N']))]])
        valA = np.append(valA,float(v_['VAL']))
        [ig,ig_,_] = s2mpj_ii('INT',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'INT')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['U'+str(int(v_['0']))]])
        valA = np.append(valA,float(0.5))
        [ig,ig_,_] = s2mpj_ii('CON'+str(int(v_['0'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'CON'+str(int(v_['0'])))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['U'+str(int(v_['0']))]])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['U'+str(int(v_['1']))]])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('CON'+str(int(v_['0'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'CON'+str(int(v_['0'])))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['W'+str(int(v_['0']))]])
        valA = np.append(valA,float(v_['-1/N**2']))
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            v_['I+1'] = 1+I
            v_['I-1'] = -1+I
            [ig,ig_,_] = s2mpj_ii('CON'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'CON'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['U'+str(I)]])
            valA = np.append(valA,float(2.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['U'+str(int(v_['I+1']))]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['U'+str(int(v_['I-1']))]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['W'+str(I)]])
            valA = np.append(valA,float(v_['-1/N**2']))
            [ig,ig_,_] = s2mpj_ii('INT',ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'INT')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['U'+str(I)]])
            valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('CON'+str(int(v_['N'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'CON'+str(int(v_['N'])))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['U'+str(int(v_['N']))]])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['U'+str(int(v_['N-1']))]])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('CON'+str(int(v_['N'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'CON'+str(int(v_['N'])))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['W'+str(int(v_['N']))]])
        valA = np.append(valA,float(v_['-1/N**2']))
        [ig,ig_,_] = s2mpj_ii('INT',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'INT')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['U'+str(int(v_['N']))]])
        valA = np.append(valA,float(0.5))
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
        self.gconst = arrset(self.gconst,ig_['INT'],float(v_['INT']))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        for I in range(int(v_['0']),int(v_['N'])+1):
            self.xlower[ix_['U'+str(I)]] = -1.0
            self.xupper[ix_['U'+str(I)]] = 1.0
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(0.0))
        for I in range(int(v_['0']),int(v_['N'])+1):
            self.x0[ix_['U'+str(I)]] = float(v_['V'+str(I)])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'Z')
        [it,iet_,_] = s2mpj_ii( 'ePROD', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['0']),int(v_['N'])+1):
            ename = 'C'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'ePROD')
            ielftype = arrset(ielftype,ie,iet_["ePROD"])
            vname = 'U'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'W'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
            posev = np.where(elftv[ielftype[ie]]=='Y')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        for I in range(int(v_['0']),int(v_['N-1'])+1):
            v_['I+1'] = 1+I
            ename = 'D'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSQ')
            ielftype = arrset(ielftype,ie,iet_["eSQ"])
            vname = 'U'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
            posev = np.where(elftv[ielftype[ie]]=='Z')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'O'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'ePROD')
            ielftype = arrset(ielftype,ie,iet_["ePROD"])
            vname = 'U'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'U'+str(int(v_['I+1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
            posev = np.where(elftv[ielftype[ie]]=='Y')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D'+str(int(v_['N']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQ')
        ielftype = arrset(ielftype,ie,iet_["eSQ"])
        ename = 'D'+str(int(v_['N']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'U'+str(int(v_['N']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='Z')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        ig = ig_['OBJ']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O'+str(int(v_['0']))])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['-TWO']))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D'+str(int(v_['0']))])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['ONE']))
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            ig = ig_['OBJ']
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['O'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['-TWO']))
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['D'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['TWO']))
        ig = ig_['OBJ']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D'+str(int(v_['N']))])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['ONE']))
        for I in range(int(v_['0']),int(v_['N'])+1):
            ig = ig_['OBJ']
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['C'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['1/N**2']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# XL SOLUTION            -4.93517D+02   $ N = 10 
# XL SOLUTION            -5.30009D+03   $ N = 100
# XL SOLUTION            -5.33674D+04   $ N = 1000
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
        self.pbclass   = "C-CQLR2-MN-V-V"
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

    @staticmethod
    def ePROD(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]
            g_[1] = EV_[0]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = 1.0
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

