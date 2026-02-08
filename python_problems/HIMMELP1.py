from s2mpjlib import *
class  HIMMELP1(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HIMMELP1
#    *********
# 
#    A nonlinear problem with inequality constraints, attributed to Himmelblau
#    by B.N. Pshenichnyj (case 0: only bounds)
# 
#    Source: 
#    B.N. Pshenichnyj
#    "The Linearization Method for Constrained Optimization",
#    Springer Verlag, SCM Series 22, Heidelberg, 1994
# 
#    SIF input: Ph. Toint, December 1994.
# 
#    classification = "C-COBR2-AN-2-0"
# 
#    Problem data
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HIMMELP1'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['B1'] = 0.1963666677
        v_['B1'] = 75.0+v_['B1']
        v_['B2'] = -.8112755343
        v_['B2'] = -3.0+v_['B2']
        v_['B6'] = -.8306567613
        v_['B6'] = -6.0+v_['B6']
        v_['-B2'] = -1.0*v_['B2']
        v_['-B6'] = -1.0*v_['B6']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [iv,ix_,_] = s2mpj_ii('X1',ix_)
        self.xnames=arrset(self.xnames,iv,'X1')
        [iv,ix_,_] = s2mpj_ii('X2',ix_)
        self.xnames=arrset(self.xnames,iv,'X2')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X1']])
        valA = np.append(valA,float(v_['-B2']))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2']])
        valA = np.append(valA,float(v_['-B6']))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        selfnob      = ngrp
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        self.gconst = arrset(self.gconst,ig_['OBJ'],float(v_['B1']))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        self.xlower[ix_['X1']] = 0.0
        self.xupper[ix_['X1']] = 95.0
        self.xlower[ix_['X2']] = 0.0
        self.xupper[ix_['X2']] = 75.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.x0[ix_['X1']] = float(95.0)
        self.x0[ix_['X2']] = float(10.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eOBNL', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        ename = 'OB'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eOBNL')
        ielftype = arrset(ielftype,ie,iet_["eOBNL"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
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
        self.grelt = loaset(self.grelt,ig,posel,ie_['OB'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN                -62.053869846
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass   = "C-COBR2-AN-2-0"
        self.objderlvl = 2

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eOBNL(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        B3 = .1269366345
        B4 = -0.20567665
        B4 = 0.01*B4
        B5 = 0.103450e-4
        B7 = .0302344793
        B8 = -0.12813448
        B8 = 0.01*B8
        B9 = 0.352599e-4
        B10 = -0.2266e-6
        B11 = 0.2564581253
        B12 = -.003460403
        B13 = 0.135139e-4
        B14 = -.1064434908
        B14 = B14-28.0
        B15 = -0.52375e-5
        B16 = -0.63e-8
        B17 = 0.7e-9
        B18 = 0.3405462
        B18 = 0.001*B18
        B19 = -0.16638e-5
        B20 = -2.86731123
        B20 = B20-0.92e-8
        A = B7*EV_[0,0]+B8*EV_[0,0]**2+B9*EV_[0,0]**3+B10*EV_[0,0]**4
        DADX = B7+2.0*B8*EV_[0,0]+3.0*B9*EV_[0,0]**2+4.0*B10*EV_[0,0]**3
        D2ADXX = 2.0*B8+6.0*B9*EV_[0,0]+12.0*B10*EV_[0,0]**2
        B = B18*EV_[0,0]+B15*EV_[0,0]**2+B16*EV_[0,0]**3
        DBDX = B18+2.0*B15*EV_[0,0]+3.0*B16*EV_[0,0]**2
        D2BDXX = 2.0*B15+6.0*B16*EV_[0,0]
        C = B3*EV_[0,0]**2+B4*EV_[0,0]**3+B5*EV_[0,0]**4
        DCDX = 2.0*B3*EV_[0,0]+3.0*B4*EV_[0,0]**2+4.0*B5*EV_[0,0]**3
        D2CDXX = 2.0*B3+6.0*B4*EV_[0,0]+12.0*B5*EV_[0,0]**2
        F = B11*EV_[1,0]**2+B12*EV_[1,0]**3+B13*EV_[1,0]**4
        DFDY = 2.0*B11*EV_[1,0]+3.0*B12*EV_[1,0]**2+4.0*B13*EV_[1,0]**3
        D2FDYY = 2.0*B11+6.0*B12*EV_[1,0]+12.0*B13*EV_[1,0]**2
        G = B17*EV_[0,0]**3+B19*EV_[0,0]
        DGDX = B19+3.0*B17*EV_[0,0]**2
        D2GDXX = 6.0*B17*EV_[0,0]
        E = np.exp(0.0005*EV_[0,0]*EV_[1,0])
        DEDX = 0.0005*EV_[1,0]*E
        DEDY = 0.0005*EV_[0,0]*E
        D2EDXX = 0.0005*EV_[1,0]*DEDX
        D2EDXY = 0.0005*(EV_[1,0]*DEDY+E)
        D2EDYY = 0.0005*EV_[0,0]*DEDY
        f_   = C+EV_[1,0]*A+F+B14/(1.0+EV_[1,0])+B*EV_[1,0]**2+G*EV_[1,0]**3+B20*E
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = DCDX+EV_[1,0]*DADX+DBDX*EV_[1,0]**2+DGDX*EV_[1,0]**3+B20*DEDX
            g_[1] = (A+DFDY-B14/(1.0+EV_[1,0])**2+2.0*B*EV_[1,0]+3.0*G*EV_[1,0]**2+
                 B20*DEDY)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = (D2CDXX+EV_[1,0]*D2ADXX+D2BDXX*EV_[1,0]**2+D2GDXX*EV_[1,0]**3+
                     B20*D2EDXX)
                H_[0,1] = DADX+2.0*EV_[1,0]*DBDX+3.0*DGDX*EV_[1,0]**2+B20*D2EDXY
                H_[1,0] = H_[0,1]
                H_[1,1] = D2FDYY+2.0*B14/(1.0+EV_[1,0])**3+2.0*B+6.0*G*EV_[1,0]+B20*D2EDYY
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

