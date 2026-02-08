from s2mpjlib import *
class  HS54(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS54
#    *********
# 
#    Source: problem 54, incorrectly stated in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    Betts problem 11.7, JOTA 21, 1977, pp.137-174.
#    SIF input: A.R. Conn, April 1990 and Nick Gould, October 1990
# 
#    classification = "C-COLR2-AN-6-1"
# 
#    some useful parameters, including N, the number of variables.
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS54'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 6
        v_['1'] = 1
        v_['6'] = 6
        v_['RHO'] = 2.0e-1
        v_['RHOSQR'] = v_['RHO']*v_['RHO']
        v_['1-RHOSQR'] = 1.0-v_['RHOSQR']
        v_['FACTOR'] = 1.0/v_['1-RHOSQR']
        v_['MU1'] = 1.0e+4
        v_['MU2'] = 1.0e+0
        v_['MU3'] = 2.0e+6
        v_['MU4'] = 1.0e+1
        v_['MU5'] = 1.0e-3
        v_['MU6'] = 1.0e+8
        v_['SIGMA1'] = 8.0e+3
        v_['SIGMA2'] = 1.0e+0
        v_['SIGMA3'] = 7.0e+6
        v_['SIGMA4'] = 5.0e+1
        v_['SIGMA5'] = 5.0e-2
        v_['SIGMA6'] = 5.0e+8
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
        [ig,ig_,_] = s2mpj_ii('CON1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'CON1')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X1']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2']])
        valA = np.append(valA,float(4.0e+3))
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
        v_['0.2SI1'] = 2.0e-1*v_['SIGMA1']
        v_['2000SI2'] = 2.0e+3*v_['SIGMA2']
        v_['4000MU2'] = 4.0e+3*v_['MU2']
        v_['RHS'] = v_['MU1']+v_['4000MU2']
        v_['RHS'] = v_['RHS']+v_['0.2SI1']
        v_['RHS'] = v_['RHS']+v_['2000SI2']
        self.gconst = arrset(self.gconst,ig_['CON1'],float(v_['RHS']))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        self.xupper[ix_['X1']] = 2.0e+4
        self.xlower[ix_['X2']] = - 1.0e+1
        self.xupper[ix_['X2']] = 1.0e+1
        self.xupper[ix_['X3']] = 1.0e+7
        self.xupper[ix_['X4']] = 2.0e+1
        self.xlower[ix_['X5']] = - 1.0e+0
        self.xupper[ix_['X5']] = 1.0e+0
        self.xupper[ix_['X6']] = 2.0e+8
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        if('X1' in ix_):
            self.x0[ix_['X1']] = float(6.0e+3)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X1']),float(6.0e+3)))
        if('X2' in ix_):
            self.x0[ix_['X2']] = float(1.5e+0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X2']),float(1.5e+0)))
        if('X3' in ix_):
            self.x0[ix_['X3']] = float(4.0e+6)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X3']),float(4.0e+6)))
        if('X4' in ix_):
            self.x0[ix_['X4']] = float(2.0e+0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X4']),float(2.0e+0)))
        if('X5' in ix_):
            self.x0[ix_['X5']] = float(3.0e-3)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X5']),float(3.0e-3)))
        if('X6' in ix_):
            self.x0[ix_['X6']] = float(5.0e+7)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X6']),float(5.0e+7)))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQR', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftp = []
        elftp = loaset(elftp,it,0,'MU')
        elftp = loaset(elftp,it,1,'SIGMA')
        [it,iet_,_] = s2mpj_ii( 'ePROD', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftp = loaset(elftp,it,0,'MU1')
        elftp = loaset(elftp,it,1,'MU2')
        elftp = loaset(elftp,it,2,'SIGMA1')
        elftp = loaset(elftp,it,3,'SIGMA2')
        elftp = loaset(elftp,it,4,'RHO')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['1']),int(v_['6'])+1):
            ename = 'E'+str(I)
            [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
            if newelt:
                self.elftype = arrset(self.elftype,ie,'eSQR')
                ielftype = arrset(ielftype,ie,iet_['eSQR'])
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='MU')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['MU'+str(I)]))
            posep = np.where(elftp[ielftype[ie]]=='SIGMA')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['SIGMA'+str(I)]))
        ename = 'F1'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePROD')
        ielftype = arrset(ielftype,ie,iet_["ePROD"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='RHO')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(v_['RHO']))
        posep = np.where(elftp[ielftype[ie]]=='MU1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(v_['MU1']))
        posep = np.where(elftp[ielftype[ie]]=='MU2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(v_['MU2']))
        posep = np.where(elftp[ielftype[ie]]=='SIGMA1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(v_['SIGMA1']))
        posep = np.where(elftp[ielftype[ie]]=='SIGMA2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(v_['SIGMA2']))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gNORMAL',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        ig = ig_['OBJ']
        self.grftype = arrset(self.grftype,ig,'gNORMAL')
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E1'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['FACTOR']))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E2'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['FACTOR']))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E3'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0e+0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E4'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0e+0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E5'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0e+0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E6'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0e+0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F1'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['FACTOR']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               0.90807482
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
        self.pbclass   = "C-COLR2-AN-6-1"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSQR(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        V1MP = (EV_[0,0]-self.elpar[iel_][0])/self.elpar[iel_][1]
        f_   = V1MP**2
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*V1MP/self.elpar[iel_][1]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0/self.elpar[iel_][1]**2
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
        TERM1 = (EV_[0,0]-self.elpar[iel_][0])/self.elpar[iel_][2]
        TERM2 = (EV_[1,0]-self.elpar[iel_][1])/self.elpar[iel_][3]
        RHO2 = self.elpar[iel_][4]+self.elpar[iel_][4]
        f_   = RHO2*TERM1*TERM2
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = RHO2*TERM2/self.elpar[iel_][2]
            g_[1] = RHO2*TERM1/self.elpar[iel_][3]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = RHO2/(self.elpar[iel_][2]*self.elpar[iel_][3])
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gNORMAL(self,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        EXPHV = np.exp(-0.5*GVAR_)
        f_= -EXPHV
        if nargout>1:
            g_ = 5.0e-1*EXPHV
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = -2.5e-1*EXPHV
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

