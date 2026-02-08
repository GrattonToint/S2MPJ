from s2mpjlib import *
class  HS68(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS68
#    *********
# 
#    This is a cost optimal inspection plan.
# 
#    Source: problem 68 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: Nick Gould, August 1991.
#               Python coding: Cunxin Huang, 2025.
# 
#    classification = "C-COOR2-MN-4-2"
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS68'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 4
        v_['A'] = 0.0001
        v_['B'] = 1.0
        v_['D'] = 1.0
        v_['NN'] = 24.0
        v_['1'] = 1
        v_['AN'] = v_['A']*v_['NN']
        v_['ROOTN'] = np.sqrt(v_['NN'])
        v_['DROOTN'] = v_['D']*v_['ROOTN']
        v_['-DROOTN'] = -v_['DROOTN']
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
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('C1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C1')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X3']])
        valA = np.append(valA,float(1.0e+0))
        [ig,ig_,_] = s2mpj_ii('C2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C2')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X4']])
        valA = np.append(valA,float(1.0e+0))
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
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        self.xlower[ix_['X1']] = 0.0001
        self.xupper[ix_['X1']] = 100.0
        self.xupper[ix_['X2']] = 100.0
        self.xupper[ix_['X3']] = 2.0
        self.xupper[ix_['X4']] = 2.0
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(1.0))
        self.y0 = np.full((self.m,1),float(1.0))
        pass
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eRECIP', iet_)
        elftv = loaset(elftv,it,0,'X1')
        [it,iet_,_] = s2mpj_ii( 'eNASTYEXP', iet_)
        elftv = loaset(elftv,it,0,'X1')
        elftv = loaset(elftv,it,1,'X3')
        elftv = loaset(elftv,it,2,'X4')
        elftp = []
        elftp = loaset(elftp,it,0,'B')
        [it,iet_,_] = s2mpj_ii( 'ePHI', iet_)
        elftv = loaset(elftv,it,0,'X2')
        elftp = loaset(elftp,it,0,'P')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        ename = 'OE1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eRECIP')
        ielftype = arrset(ielftype,ie,iet_["eRECIP"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(1.0))
        posev = np.where(elftv[ielftype[ie]]=='X1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'OE2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eNASTYEXP')
        ielftype = arrset(ielftype,ie,iet_["eNASTYEXP"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(1.0))
        posev = np.where(elftv[ielftype[ie]]=='X1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(1.0))
        posev = np.where(elftv[ielftype[ie]]=='X3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(1.0))
        posev = np.where(elftv[ielftype[ie]]=='X4')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='B')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(v_['B']))
        ename = 'C1E1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePHI')
        ielftype = arrset(ielftype,ie,iet_["ePHI"])
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(1.0))
        posev = np.where(elftv[ielftype[ie]]=='X2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='P')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(0.0e+0))
        ename = 'C2E1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePHI')
        ielftype = arrset(ielftype,ie,iet_["ePHI"])
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(1.0))
        posev = np.where(elftv[ielftype[ie]]=='X2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='P')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(v_['DROOTN']))
        ename = 'C2E2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePHI')
        ielftype = arrset(ielftype,ie,iet_["ePHI"])
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(1.0))
        posev = np.where(elftv[ielftype[ie]]=='X2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='P')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(v_['-DROOTN']))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        ig = ig_['OBJ']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['OE1'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['AN']))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['OE2'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0e+0))
        ig = ig_['C1']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['C1E1'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-2.0e+0))
        ig = ig_['C2']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['C2E1'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0e+0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['C2E2'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0e+0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# XL SOLUTION            -9.20389D-01
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
        self.pbclass   = "C-COOR2-MN-4-2"
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
        self.efpar = arrset( self.efpar,0,3.9894228040143270e-01)
        return pbm

    @staticmethod
    def eRECIP(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        x    = EV_[0].item()
        iel_ = args[1]
        f_   = 1.0e+0/x
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -1.0e+0/x**2
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0e+0/x**3
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eNASTYEXP(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        x    = EV_[0].item()
        y    = EV_[1].item()
        z    = EV_[2].item()
        iel_ = args[1]
        E = np.exp(x)
        R = self.elpar[iel_][0]*(E-1.0e+0)-y
        S = E-1.0e+0+z
        F1 = -(-(R*z/S))/(x*x)
        F2 = -R/(S*x)
        F3 = -z/(S*x)
        F4 = (z*R)/(x*S*S)
        F11 = 2.0e+0*(-(R*z/S))/(x*x*x)
        F12 = E*z*((R/S)-self.elpar[iel_][0])/(x*S)
        F11 = F11+F12
        F12 = R/(S*x*x)
        F13 = z/(S*x*x)
        F14 = R/(S*S*x)
        F15 = -1.0e+0/(S*x)
        F16 = -z*R/(S*S*x*x)
        F17 = z/(S*S*x)
        F18 = -2.0e+0*z*R/(x*S*S*S)
        f_   = z*R/(S*x)
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -F1-(self.elpar[iel_][0]*E*F3)-(E*F4)
            g_[1] = F3
            g_[2] = -F2-F4
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = (-F11-(2.0e+0*self.elpar[iel_][0]*E*F13)-(2.0e+0*E*F16)-
                     (2.0e+0*self.elpar[iel_][0]*E*E*F17)-(E*E*F18))
                H_[1,0] = F13+(E*F17)
                H_[0,1] = H_[1,0]
                H_[2,0] = (-(E*F14)-(self.elpar[iel_][0]*E*F15)-F16-(self.elpar[iel_][0]*E*F17)-
                     (E*F18)-F12)
                H_[0,2] = H_[2,0]
                H_[2,1] = F15+F17
                H_[1,2] = H_[2,1]
                H_[2,2] = -(2.0e+0*F14)-F18
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ePHI(self, nargout,*args):

        import numpy as np
        from scipy import special
        x2    = args[0]
        iel_  = args[1]
        mx2pp = -x2.item()+self.elpar[iel_][0]
        E     = np.exp(-5.0e-1*mx2pp**2)
        arg   =  7.071067811865475e-1 * np.abs( mx2pp )
        if mx2pp >= 0 :
            f_ =  0.5 + 0.5 * special.erf( arg )
        else:
            f_ =  0.5 - 0.5 * special.erf( arg )
        if nargout>1:
            dim  = len(x2)
            g_    = np.zeros(dim)
            g_[0] = -self.efpar[0]*E
            if nargout>2:
                H_      = np.zeros((1,1))
                H_[0,0] = -self.efpar[0]*E*mx2pp
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
