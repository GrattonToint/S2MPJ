from s2mpjlib import *
class  LUKSAN12(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : LUKSAN12
#    *********
# 
#    Problem 12 (chained and modified HS47) in the paper
# 
#      L. Luksan
#      Hybrid methods in large sparse nonlinear least squares
#      J. Optimization Theory & Applications 89(3) 575-595 (1996)
# 
#    SIF input: Nick Gould, June 2017.
# 
#    classification = "C-CNOR2-AN-V-V"
# 
#   seed for dimensions
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 17 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'LUKSAN12'

    def __init__(self, *args): 
        import numpy as np
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['S'] = 32
        v_['N'] = 3*v_['S']
        v_['N'] = 2+v_['N']
        v_['M'] = 6*v_['S']
        v_['1'] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.A       = lil_matrix((1000000,1000000))
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames      = np.array([])
        self.cnames = np.array([])
        gtype       = np.array([])
        v_['I'] = 1
        v_['K'] = 1
        for J in range(int(v_['1']),int(v_['S'])+1):
            v_['K+1'] = 1+v_['K']
            v_['K+2'] = 2+v_['K']
            v_['K+3'] = 3+v_['K']
            v_['K+4'] = 4+v_['K']
            v_['K+5'] = 5+v_['K']
            v_['I+1'] = 1+v_['I']
            v_['I+2'] = 2+v_['I']
            [ig,ig_,_] = s2mpj_ii('E'+str(int(v_['K'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'E'+str(int(v_['K'])))
            iv = ix_['X'+str(int(v_['I+1']))]
            self.A[ig,iv] = float(-10.0e0)+self.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('E'+str(int(v_['K+1'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'E'+str(int(v_['K+1'])))
            iv = ix_['X'+str(int(v_['I+2']))]
            self.A[ig,iv] = float(1.0e0)+self.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('E'+str(int(v_['K+2'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'E'+str(int(v_['K+2'])))
            [ig,ig_,_] = s2mpj_ii('E'+str(int(v_['K+3'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'E'+str(int(v_['K+3'])))
            [ig,ig_,_] = s2mpj_ii('E'+str(int(v_['K+4'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'E'+str(int(v_['K+4'])))
            [ig,ig_,_] = s2mpj_ii('E'+str(int(v_['K+5'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'E'+str(int(v_['K+5'])))
            iv = ix_['X'+str(int(v_['I+1']))]
            self.A[ig,iv] = float(1.0e0)+self.A[ig,iv]
            v_['I'] = 3+v_['I']
            v_['K'] = 6+v_['K']
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
        self.cnames= cnames[self.congrps]
        self.nob = ngrp-self.m
        self.objgrps = np.where(gtype=='<>')[0]
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        v_['K'] = 1
        for J in range(int(v_['1']),int(v_['S'])+1):
            v_['K+1'] = 1+v_['K']
            v_['K+4'] = 4+v_['K']
            v_['K+5'] = 5+v_['K']
            self.gconst = arrset(self.gconst,ig_['E'+str(int(v_['K+1']))],float(1.0))
            self.gconst = arrset(self.gconst,ig_['E'+str(int(v_['K+4']))],float(10.0))
            self.gconst = arrset(self.gconst,ig_['E'+str(int(v_['K+5']))],float(20.0))
            v_['K'] = 6+v_['K']
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        self.xlower = np.zeros((self.n,1))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        for I in range(int(v_['1']),int(v_['N'])+1):
            self.x0[ix_['X'+str(I)]] = float(-1.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eE1', iet_)
        elftv = loaset(elftv,it,0,'X0')
        [it,iet_,_] = s2mpj_ii( 'eE3', iet_)
        elftv = loaset(elftv,it,0,'X3')
        [it,iet_,_] = s2mpj_ii( 'eE4', iet_)
        elftv = loaset(elftv,it,0,'X4')
        [it,iet_,_] = s2mpj_ii( 'eE5', iet_)
        elftv = loaset(elftv,it,0,'X0')
        elftv = loaset(elftv,it,1,'X3')
        [it,iet_,_] = s2mpj_ii( 'eF5', iet_)
        elftv = loaset(elftv,it,0,'X3')
        elftv = loaset(elftv,it,1,'X4')
        [it,iet_,_] = s2mpj_ii( 'eE6', iet_)
        elftv = loaset(elftv,it,0,'X2')
        elftv = loaset(elftv,it,1,'X3')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        v_['I'] = 1
        v_['K'] = 1
        for J in range(int(v_['1']),int(v_['S'])+1):
            v_['K+2'] = 2+v_['K']
            v_['K+3'] = 3+v_['K']
            v_['K+4'] = 4+v_['K']
            v_['K+5'] = 5+v_['K']
            v_['I+2'] = 2+v_['I']
            v_['I+3'] = 3+v_['I']
            v_['I+4'] = 4+v_['I']
            ename = 'E'+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eE1')
            ielftype = arrset(ielftype,ie,iet_["eE1"])
            ename = 'E'+str(int(v_['K']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(int(v_['I']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='X0')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'E'+str(int(v_['K+2']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eE3')
            ielftype = arrset(ielftype,ie,iet_["eE3"])
            ename = 'E'+str(int(v_['K+2']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(int(v_['I+3']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='X3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'E'+str(int(v_['K+3']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eE4')
            ielftype = arrset(ielftype,ie,iet_["eE4"])
            ename = 'E'+str(int(v_['K+3']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(int(v_['I+4']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='X4')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'E'+str(int(v_['K+4']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eE5')
            ielftype = arrset(ielftype,ie,iet_["eE5"])
            ename = 'E'+str(int(v_['K+4']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(int(v_['I']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='X0')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'E'+str(int(v_['K+4']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(int(v_['I+3']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='X3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'F'+str(int(v_['K+4']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eF5')
            ielftype = arrset(ielftype,ie,iet_["eF5"])
            ename = 'F'+str(int(v_['K+4']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(int(v_['I+3']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='X3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'F'+str(int(v_['K+4']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(int(v_['I+4']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='X4')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'E'+str(int(v_['K+5']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eE6')
            ielftype = arrset(ielftype,ie,iet_["eE6"])
            ename = 'E'+str(int(v_['K+5']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(int(v_['I+2']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='X2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'E'+str(int(v_['K+5']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(int(v_['I+3']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='X3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            v_['I'] = 3+v_['I']
            v_['K'] = 6+v_['K']
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        v_['K'] = 1
        for J in range(int(v_['1']),int(v_['S'])+1):
            v_['K+2'] = 2+v_['K']
            v_['K+3'] = 3+v_['K']
            v_['K+4'] = 4+v_['K']
            v_['K+5'] = 5+v_['K']
            ig = ig_['E'+str(int(v_['K']))]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['K']))])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
            ig = ig_['E'+str(int(v_['K+2']))]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['K+2']))])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
            ig = ig_['E'+str(int(v_['K+3']))]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['K+3']))])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
            ig = ig_['E'+str(int(v_['K+4']))]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['K+4']))])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = posel+1
            self.grelt = loaset(self.grelt,ig,posel,ie_['F'+str(int(v_['K+4']))])
            self.grelw = loaset(self.grelw,ig,posel, 1.)
            ig = ig_['E'+str(int(v_['K+5']))]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['K+5']))])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
            v_['K'] = 6+v_['K']
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 0.0
#    Solution
# LO SOLTN                0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        self.A.resize(ngrp,self.n)
        self.A     = self.A.tocsr()
        sA1,sA2    = self.A.shape
        self.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass = "C-CNOR2-AN-V-V"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eE1(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = 10.0e0*EV_[0]**2
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 20.0e0*EV_[0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 20.0e0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eE3(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = (EV_[0]-1.0e0)**2
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0e0*(EV_[0]-1.0e0)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0e0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eE4(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = (EV_[0]-1.0e0)**3
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 3.0e0*(EV_[0]-1.0e0)**2
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 6.0e0*(EV_[0]-1.0e0)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eE5(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[1]*EV_[0]*EV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0e0*EV_[1]*EV_[0]
            g_[1] = EV_[0]*EV_[0]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 2.0e0*EV_[1]
                H_[0,1] = 2.0e0*EV_[0]
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eF5(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((1,2))
        IV_ = np.zeros(1)
        U_[0,0] = U_[0,0]+1
        U_[0,1] = U_[0,1]-1
        IV_[0] = U_[0:1,:].dot(EV_)
        f_   = np.sin(IV_[0])
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = np.cos(IV_[0])
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = -np.sin(IV_[0])
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eE6(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = (EV_[0]**4)*(EV_[1]**2)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 4.0e0*(EV_[0]**3)*(EV_[1]**2)
            g_[1] = 2.0e0*(EV_[0]**4)*EV_[1]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 12.0e0*(EV_[0]**2)*(EV_[1]**2)
                H_[0,1] = 8.0e0*(EV_[0]**3)*EV_[1]
                H_[1,0] = H_[0,1]
                H_[1,1] = 2.0e0*(EV_[0]**4)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

