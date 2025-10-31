from s2mpjlib import *
class  MINSURF(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : MINSURF
#    *********
#    Variable dimension full rank linear problem
#    A version of the minimum surface problem
#    on the unit square with simple boundary conditions.
# 
#    SIF input: Ph. Toint, Jan 1991.
# 
#    classification = "C-COXR2-MY-64-0"
# 
#    Discretization parameter
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'MINSURF'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['P'] = 7
        v_['1'] = 1
        v_['P+1'] = 1+v_['P']
        v_['RP'] = float(v_['P'])
        v_['RPSQ'] = v_['RP']*v_['RP']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for i in range(int(v_['1']),int(v_['P+1'])+1):
            for j in range(int(v_['1']),int(v_['P+1'])+1):
                [iv,ix_,_] = s2mpj_ii('X'+str(i)+','+str(j),ix_)
                self.xnames=arrset(self.xnames,iv,'X'+str(i)+','+str(j))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for i in range(int(v_['1']),int(v_['P'])+1):
            for j in range(int(v_['1']),int(v_['P'])+1):
                [ig,ig_,_] = s2mpj_ii('S'+str(i)+','+str(j),ig_)
                gtype = arrset(gtype,ig,'<>')
                self.gscale = arrset(self.gscale,ig,float(v_['RPSQ']))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        for i in range(int(v_['1']),int(v_['P'])+1):
            for j in range(int(v_['1']),int(v_['P'])+1):
                self.gconst = arrset(self.gconst,ig_['S'+str(i)+','+str(j)],float(-1.0))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        v_['2'] = 2
        for i in range(int(v_['2']),int(v_['P'])+1):
            for j in range(int(v_['2']),int(v_['P'])+1):
                self.xlower[ix_['X'+str(i)+','+str(j)]] = -float('Inf')
                self.xupper[ix_['X'+str(i)+','+str(j)]] = +float('Inf')
        for i in range(int(v_['1']),int(v_['P+1'])+1):
            self.xlower[ix_['X'+str(int(v_['1']))+','+str(i)]] = 1.0
            self.xupper[ix_['X'+str(int(v_['1']))+','+str(i)]] = 1.0
            self.xlower[ix_['X'+str(int(v_['P+1']))+','+str(i)]] = 1.0
            self.xupper[ix_['X'+str(int(v_['P+1']))+','+str(i)]] = 1.0
            self.xlower[ix_['X'+str(i)+','+str(int(v_['1']))]] = 1.0
            self.xupper[ix_['X'+str(i)+','+str(int(v_['1']))]] = 1.0
            self.xlower[ix_['X'+str(i)+','+str(int(v_['P+1']))]] = 1.0
            self.xupper[ix_['X'+str(i)+','+str(int(v_['P+1']))]] = 1.0
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eISQ', iet_)
        elftv = loaset(elftv,it,0,'V')
        elftv = loaset(elftv,it,1,'W')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for i in range(int(v_['1']),int(v_['P'])+1):
            v_['i+1'] = 1+i
            for j in range(int(v_['1']),int(v_['P'])+1):
                v_['j+1'] = 1+j
                ename = 'A'+str(i)+','+str(j)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eISQ')
                ielftype = arrset(ielftype,ie,iet_["eISQ"])
                self.x0 = np.zeros((self.n,1))
                vname = 'X'+str(i)+','+str(j)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                posev = np.where(elftv[ielftype[ie]]=='V')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'X'+str(int(v_['i+1']))+','+str(int(v_['j+1']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                posev = np.where(elftv[ielftype[ie]]=='W')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                ename = 'B'+str(i)+','+str(j)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eISQ')
                ielftype = arrset(ielftype,ie,iet_["eISQ"])
                vname = 'X'+str(i)+','+str(int(v_['j+1']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                posev = np.where(elftv[ielftype[ie]]=='V')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'X'+str(int(v_['i+1']))+','+str(j)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                posev = np.where(elftv[ielftype[ie]]=='W')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gSQROOT',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        v_['WEIGHT'] = 0.5*v_['RPSQ']
        for i in range(int(v_['1']),int(v_['P'])+1):
            for j in range(int(v_['1']),int(v_['P'])+1):
                ig = ig_['S'+str(i)+','+str(j)]
                self.grftype = arrset(self.grftype,ig,'gSQROOT')
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['A'+str(i)+','+str(j)])
                self.grelw = loaset(self.grelw,ig,posel,float(v_['WEIGHT']))
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['B'+str(i)+','+str(j)])
                self.grelw = loaset(self.grelw,ig,posel,float(v_['WEIGHT']))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass   = "C-COXR2-MY-64-0"
        self.x0        = np.zeros((self.n,1))
        self.objderlvl = 2


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eISQ(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((1,2))
        IV_ = np.zeros(1)
        U_[0,0] = U_[0,0]+1
        U_[0,1] = U_[0,1]-1
        IV_[0] = U_[0:1,:].dot(EV_)
        f_   = IV_[0]*IV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = IV_[0]+IV_[0]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gSQROOT(self,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        SQRAL = np.sqrt(GVAR_)
        f_= SQRAL
        if nargout>1:
            g_ = 0.5/SQRAL
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = -0.25/(GVAR_*SQRAL)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

