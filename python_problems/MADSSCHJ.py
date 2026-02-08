from s2mpjlib import *
class  MADSSCHJ(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : MADSSCHJ
#    *********
# 
#    A nonlinear minmax problem with variable dimension.
#    The Jacobian of the constraints is dense.
# 
#    Source:
#    K. Madsen and H. Schjaer-Jacobsen,
#    "Linearly Constrained Minmax Optimization",
#    Mathematical Programming 14, pp. 208-223, 1978.
# 
#    SIF input: Ph. Toint, August 1993.
# 
#    classification = "C-CLQR2-AN-V-V"
# 
#    N is the number of variables - 1, and must be even and at least 4.
#    The number of inequality constraints is 2*N - 2.
# 
#           Alternative values for the SIF file parameters:
# IE N                   4              $-PARAMETER  n=  5, m=  6
# IE N                   10             $-PARAMETER  n= 11, m= 18  original value
# IE N                   20             $-PARAMETER  n= 21, m= 38
# IE N                   30             $-PARAMETER  n= 31, m= 58
# IE N                   40             $-PARAMETER  n= 41, m= 78
# IE N                   50             $-PARAMETER  n= 51, m= 98
# IE N                   60             $-PARAMETER  n= 61, m=118
# IE N                   70             $-PARAMETER  n= 71, m=138
# IE N                   80             $-PARAMETER  n= 81, m=158
# IE N                   90             $-PARAMETER  n= 91, m=178
# IE N                   100            $-PARAMETER  n=101, m=198
# IE N                   200            $-PARAMETER  n=201, m=398
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'MADSSCHJ'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(4);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
        v_['1'] = 1
        v_['2'] = 2
        v_['3'] = 3
        v_['4'] = 4
        v_['N-1'] = -1+v_['N']
        v_['2N'] = v_['N']+v_['N']
        v_['M'] = -2+v_['2N']
        v_['M-1'] = -1+v_['M']
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
        [iv,ix_,_] = s2mpj_ii('Z',ix_)
        self.xnames=arrset(self.xnames,iv,'Z')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Z']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('C1',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C1')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Z']])
        valA = np.append(valA,float(1.0))
        for I in range(int(v_['2']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('C1',ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'C1')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(I)]])
            valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('C2',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C2')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Z']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X1']])
        valA = np.append(valA,float(-1.0))
        for I in range(int(v_['3']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('C2',ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'C2')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(I)]])
            valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('C3',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C3')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Z']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X1']])
        valA = np.append(valA,float(-1.0))
        for I in range(int(v_['3']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('C3',ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'C3')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(I)]])
            valA = np.append(valA,float(-1.0))
        for K in range(int(v_['4']),int(v_['M-1'])+1,int(v_['2'])):
            v_['K+1'] = 1+K
            v_['K+2'] = 2+K
            v_['J'] = int(np.fix(v_['K+2']/v_['2']))
            v_['J-1'] = -1+v_['J']
            v_['J+1'] = 1+v_['J']
            [ig,ig_,_] = s2mpj_ii('C'+str(K),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'C'+str(K))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['Z']])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('C'+str(int(v_['K+1'])),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'C'+str(int(v_['K+1'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['Z']])
            valA = np.append(valA,float(1.0))
            for I in range(int(v_['1']),int(v_['J-1'])+1):
                [ig,ig_,_] = s2mpj_ii('C'+str(K),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'C'+str(K))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(I)]])
                valA = np.append(valA,float(-1.0))
                [ig,ig_,_] = s2mpj_ii('C'+str(int(v_['K+1'])),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'C'+str(int(v_['K+1'])))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(I)]])
                valA = np.append(valA,float(-1.0))
            for I in range(int(v_['J+1']),int(v_['N'])+1):
                [ig,ig_,_] = s2mpj_ii('C'+str(K),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'C'+str(K))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(I)]])
                valA = np.append(valA,float(-1.0))
                [ig,ig_,_] = s2mpj_ii('C'+str(int(v_['K+1'])),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'C'+str(int(v_['K+1'])))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(I)]])
                valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('C'+str(int(v_['M'])),ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C'+str(int(v_['M'])))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Z']])
        valA = np.append(valA,float(1.0))
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            [ig,ig_,_] = s2mpj_ii('C'+str(int(v_['M'])),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'C'+str(int(v_['M'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(I)]])
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
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        for K in range(int(v_['1']),int(v_['M'])+1):
            self.gconst = arrset(self.gconst,ig_['C'+str(K)],float(-1.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(10.0))
        self.x0[ix_['Z']] = float(0.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['1']),int(v_['N'])+1):
            ename = 'XSQ'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSQ')
            ielftype = arrset(ielftype,ie,iet_["eSQ"])
            vname = 'X'+str(I)
            [iv,ix_]  = (
                  s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(10.0)))
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        ig = ig_['C'+str(int(v_['1']))]
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['XSQ'+str(int(v_['1']))])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        ig = ig_['C'+str(int(v_['2']))]
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['XSQ'+str(int(v_['2']))])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        ig = ig_['C'+str(int(v_['3']))]
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['XSQ'+str(int(v_['2']))])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-2.0))
        for K in range(int(v_['4']),int(v_['M-1'])+1,int(v_['2'])):
            v_['K+1'] = 1+K
            v_['K+2'] = 2+K
            v_['J'] = int(np.fix(v_['K+2']/v_['2']))
            ig = ig_['C'+str(K)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['XSQ'+str(int(v_['J']))])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
            ig = ig_['C'+str(int(v_['K+1']))]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['XSQ'+str(int(v_['J']))])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(-2.0))
        ig = ig_['C'+str(int(v_['M']))]
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['XSQ'+str(int(v_['N']))])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN(4)             -2.6121094144
# LO SOLTN(10)            -12.814452425
# LO SOLTN(20)            -49.869888156
# LO SOLTN(30)            -111.93545559
# LO SOLTN(40)            -199.00371592
# LO SOLTN(50)            -311.07308068
# LO SOLTN(60)            -448.14300524
# LO SOLTN(70)            -610.21325256
# LO SOLTN(80)            -797.28370289
# LO SOLTN(90)            -1009.3542892
# LO SOLTN(100)           -1246.4249710
# LO SOLTN(200)           -4992.1339031
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle+self.neq,self.m)] = np.zeros((self.nge,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CLQR2-AN-V-V"
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

