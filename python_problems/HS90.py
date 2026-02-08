from s2mpjlib import *
class  HS90(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS90
#    *********
# 
#    A time-optimal heat conduction problem.
# 
#    Source: problem 91 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: Nick Gould, September 1991.
#      SAVEs removed December 3rd 2014
#      Python coding: Cunxin Huang, 2025.
# 
#    classification = "C-CQOR2-MN-6-1"
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS90'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 6
        v_['EPS'] = 0.01
        v_['EPSSQR'] = v_['EPS']*v_['EPS']
        v_['-EPSSQR'] = -1.0*v_['EPSSQR']
        v_['1'] = 1
        v_['2'] = 2
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
        [ig,ig_,_] = s2mpj_ii('CON',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'CON')
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
        self.gconst = arrset(self.gconst,ig_['CON'],float(v_['-EPSSQR']))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        for I in range(int(v_['1']),int(v_['N'])+1,int(v_['2'])):
            if('X'+str(I) in ix_):
                self.x0[ix_['X'+str(I)]] = float(0.5)
            else:
                self.y0  = (
                      arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X'+str(I)]),float(0.5)))
        for I in range(int(v_['2']),int(v_['N'])+1,int(v_['2'])):
            if('X'+str(I) in ix_):
                self.x0[ix_['X'+str(I)]] = float(-0.5)
            else:
                self.y0  = (
                      arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X'+str(I)]),float(-0.5)))
        pass
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQR', iet_)
        elftv = loaset(elftv,it,0,'X')
        [it,iet_,_] = s2mpj_ii( 'eH', iet_)
        elftv = loaset(elftv,it,0,'X1')
        elftv = loaset(elftv,it,1,'X2')
        elftv = loaset(elftv,it,2,'X3')
        elftv = loaset(elftv,it,3,'X4')
        elftv = loaset(elftv,it,4,'X5')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['1']),int(v_['N'])+1):
            ename = 'O'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_["eSQR"])
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'H'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eH')
        ielftype = arrset(ielftype,ie,iet_["eH"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X4')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X5')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            ig = ig_['OBJ']
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['O'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['CON']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['H'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle+self.neq,self.m)] = np.zeros((self.nge,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CQOR2-MN-5-1"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%
    
    @staticmethod
    def extfunc(self,x):
        # A translation of the Fortran code present in the SIF file.
        import numpy as np
    

        n     = len(x)
        g     = np.zeros((n,1))
        H     = np.zeros((n,n))
        A     = np.zeros((30,1))
        R     = np.zeros((30,30))
        S     = np.zeros((30,1))
        RHO   = np.zeros((30,1))
        DRHO  = np.zeros((30,n))
        D2RHO = np.zeros((30,n,n))
        P     = np.zeros((n+1,1))
        
        mu = np.array([ 8.6033358901938017e-01, 3.4256184594817283e+00, 6.4372981791719468e+00, 9.5293344053619631e+00,
                        1.2645287223856643e+01, 1.5771284874815882e+01, 1.8902409956860023e+01, 2.2036496727938566e+01, 
                        2.5172446326646664e+01, 2.8309642854452012e+01, 3.1447714637546234e+01, 3.4586424215288922e+01, 
                        3.7725612827776501e+01, 4.0865170330488070e+01, 4.4005017920830845e+01, 4.7145097736761031e+01, 
                        5.0285366337773652e+01, 5.3425790477394663e+01, 5.6566344279821521e+01, 5.9707007305335459e+01, 
                        6.2847763194454451e+01, 6.5988598698490392e+01, 6.9129502973895256e+01, 7.2270467060308960e+01, 
                        7.5411483488848148e+01, 7.8552545984242926e+01, 8.1693649235601683e+01, 8.4834788718042290e+01, 
                        8.7975960552493220e+01, 9.1117161394464745e+01 ] )

        T = 2.0 / 15.0
        for i in np.arange(0,30):
            MUI    = mu[i]
            SMUI   = np.sin(MUI)
            CMUI   = np.cos(MUI)
            AI     = 2.0*SMUI/(MUI+SMUI*CMUI)
            A[i]   = AI
            S[i]   = 2.0*AI*(CMUI-SMUI/MUI)
            AIMUI2 = AI*MUI**2
            for j in np.arange(i+1):
                if i == j:
                    R[i,i] = 0.5*(1.0+0.5*np.sin(MUI+MUI)/MUI)*AIMUI2**2
                else:
                    MUJ    = mu[j]
                    R[i,j] = 0.5*(np.sin(MUI+MUJ )/(MUI+MUJ)+np.sin(MUI-MUJ )/(MUI-MUJ))*AIMUI2*A[j,0]*MUJ**2
                    R[j,i] = R[i,j]
        
        #                                  n   2
        #  Calculate the functions p(x) = SUM x .
        #                           j     i=j  i

        for k in np.arange(n-1,-1,-1):
            P[k] = P[k+1]+(x[k].item())**2

        #  Calculate the functions rho.

        for j in np.arange(30):
            MUJ2 = mu[j]*mu[j]
            U    = np.exp(-MUJ2*P[0,0])
            for k in np.arange(n):
                 DRHO[j,k] = 2.0*U*x[k].item()
                 for l in np.arange(k,n):
                    D2RHO[j,k,l] = -4.0*MUJ2*U*x[k].item()*x[l].item()
                    if l == k:
                        D2RHO[j,k,l] = D2RHO[j,k,l]+2.0*U
            ALPHA = -2.0
            for i in np.arange(1,n):
                EU = ALPHA*np.exp(-MUJ2*P[i,0])
                U  = U+EU
                for k in np.arange(i,n):
                     DRHO[j,k] = DRHO[j,k]+2.0*EU*x[k].item()
                     for l in range(k,n):
                         D2RHO[j,k,l] = D2RHO[j,k,l]-4.0*MUJ2*EU*x[k].item()*x[l].item()
                         if l == k :
                             D2RHO[j,k,l] = D2RHO[j,k,l]+2.0*EU
                ALPHA = -ALPHA
            U      = U+0.5*ALPHA
            RHO[j] = -U/MUJ2
            
        #  Evaluate the function and derivatives.

        f = T;
        for i in np.arange(30):
            SI   = S[i,0]
            RHOI = RHO[i,0]
            f    = f+SI*RHOI
            for k in np.arange(n):
                g[k] = g[k,0]+SI*DRHO[i,k]
                for l in np.arange(k,n):
                     H[k,l] = H[k,l]+SI*D2RHO[i,k,l]
            for j in np.arange(30):
                RIJ  = R[i,j]
                RHOJ = RHO[j,0]
                f    = f+RIJ*RHOI*RHOJ
                for k in np.arange(n):
                    g[k]= g[k,0]+RIJ*(RHOI*DRHO[j,k]+RHOJ*DRHO[i,k])
                    for l in range(k,n):
                        H[k,l] = H[k,l]+RIJ*(RHOI*D2RHO[j,k,l]+RHOJ*D2RHO[i,k,l]+DRHO[i,k]*DRHO[j,l]+DRHO[j,k]*DRHO[i,l])

        #   Symmetrize the Hessian.

        for k in np.arange(n):
            for l in np.arange(k+1,n):
                H[l,k] = H[k,l]
 
        return f, g, H
 
    @staticmethod
    def eSQR(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        x    = EV_[0,0]
        f_   = x*x
        if nargout>1:
            dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = x+x
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0e+0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eH(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        f_, g_, H_ = self.extfunc(self, EV_ )
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

