from s2mpjlib import *
class  OET5(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : OET5
#    *********
# 
#    A nonlinear programming formulation of a discretization of
#    a nonlinear Chebychev problem.
# 
#    The problem is
# 
#        min  max || phi(x,w) ||, for all w in the interval I.
#         x    w
# 
#    I is discretized, and the problem solved over the
#    discrete points.
# 
#    Nonlinear programming formulation
#        min   u     s.t.  u - phi >= 0, u + phi >= 0
#        x,u
# 
#    Specific problem: I = [0.25,1]
#    phi(x,w) = sqrt(w) - x4 - (x1 w^2 + x2 w + x3 )^2
# 
#    Source: K. Oettershagen
#    "Ein superlinear knonvergenter algorithmus zur losung 
#     semi-infiniter optimierungsproblem",
#     Ph.D thesis, Bonn University, 1982
# 
#    SIF input: Nick Gould, February, 1994.
# 
#    classification = "C-CLQR2-AN-5-V"
# 
#    Discretization
# 
# IE M                   2
# IE M                   100
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 25 XI 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'OET5'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['M'] = 500
        v_['LOWER'] = 0.25
        v_['UPPER'] = 1.0
        v_['0'] = 0
        v_['1'] = 1
        v_['DIFF'] = v_['UPPER']-v_['LOWER']
        v_['RM'] = float(v_['M'])
        v_['H'] = v_['DIFF']/v_['RM']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [iv,ix_,_] = s2mpj_ii('U',ix_)
        self.xnames=arrset(self.xnames,iv,'U')
        [iv,ix_,_] = s2mpj_ii('X1',ix_)
        self.xnames=arrset(self.xnames,iv,'X1')
        [iv,ix_,_] = s2mpj_ii('X2',ix_)
        self.xnames=arrset(self.xnames,iv,'X2')
        [iv,ix_,_] = s2mpj_ii('X3',ix_)
        self.xnames=arrset(self.xnames,iv,'X3')
        [iv,ix_,_] = s2mpj_ii('X4',ix_)
        self.xnames=arrset(self.xnames,iv,'X4')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['U']])
        valA = np.append(valA,float(1.0))
        for I in range(int(v_['0']),int(v_['M'])+1):
            [ig,ig_,_] = s2mpj_ii('LO'+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'LO'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['U']])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X4']])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('UP'+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'UP'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['U']])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X4']])
            valA = np.append(valA,float(1.0))
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
        for I in range(int(v_['0']),int(v_['M'])+1):
            v_['RI'] = float(I)
            v_['W'] = v_['RI']*v_['H']
            v_['W'] = v_['W']+v_['LOWER']
            v_['SQRTW'] = np.sqrt(v_['W'])
            v_['-SQRTW'] = -1.0*v_['SQRTW']
            self.gconst = arrset(self.gconst,ig_['LO'+str(I)],float(v_['-SQRTW']))
            self.gconst = arrset(self.gconst,ig_['UP'+str(I)],float(v_['SQRTW']))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(1.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eQUAD', iet_)
        elftv = loaset(elftv,it,0,'X1')
        elftv = loaset(elftv,it,1,'X2')
        elftv = loaset(elftv,it,2,'X3')
        elftp = []
        elftp = loaset(elftp,it,0,'W')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['0']),int(v_['M'])+1):
            v_['RI'] = float(I)
            v_['W'] = v_['RI']*v_['H']
            v_['W'] = v_['W']+v_['LOWER']
            ename = 'E'+str(int(v_['1']))+','+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eQUAD')
            ielftype = arrset(ielftype,ie,iet_["eQUAD"])
            ename = 'E'+str(int(v_['1']))+','+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X1'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(1.0))
            posev = np.where(elftv[ielftype[ie]]=='X1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'E'+str(int(v_['1']))+','+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(1.0))
            posev = np.where(elftv[ielftype[ie]]=='X2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'E'+str(int(v_['1']))+','+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(1.0))
            posev = np.where(elftv[ielftype[ie]]=='X3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'E'+str(int(v_['1']))+','+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            posep = np.where(elftp[ielftype[ie]]=='W')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['W']))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['0']),int(v_['M'])+1):
            ig = ig_['LO'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt  = (
                  loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['1']))+','+str(I)]))
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
            ig = ig_['UP'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt  = (
                  loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['1']))+','+str(I)]))
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
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
        self.pbclass   = "C-CLQR2-AN-5-V"
        self.objderlvl = 2
        self.conderlvl = [2]


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eQUAD(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        W2 = self.elpar[iel_][0]*self.elpar[iel_][0]
        TERM = W2*EV_[0]+self.elpar[iel_][0]*EV_[1]+EV_[2]
        f_   = TERM*TERM
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*TERM*W2
            g_[1] = 2.0*TERM*self.elpar[iel_][0]
            g_[2] = 2.0*TERM
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = 2.0*W2*W2
                H_[1,0] = 2.0*W2*self.elpar[iel_][0]
                H_[0,1] = H_[1,0]
                H_[2,0] = 2.0*W2
                H_[0,2] = H_[2,0]
                H_[1,1] = 2.0*W2
                H_[2,1] = 2.0*self.elpar[iel_][0]
                H_[1,2] = H_[2,1]
                H_[2,2] = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

