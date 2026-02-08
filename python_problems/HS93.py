from s2mpjlib import *
class  HS93(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS93
#    *********
# 
#    A transformer design problem.
# 
#    Source: problem 93 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: Nick Gould, August 1991.
# 
#    classification = "C-COOR2-MY-6-2"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS93'

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
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C1')
        [ig,ig_,_] = s2mpj_ii('C2',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'C2')
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
        self.gconst = arrset(self.gconst,ig_['C1'],float(2.07e+0))
        self.gconst = arrset(self.gconst,ig_['C2'],float(1.0e+0))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        if('X1' in ix_):
            self.x0[ix_['X1']] = float(5.54)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X1']),float(5.54)))
        if('X2' in ix_):
            self.x0[ix_['X2']] = float(4.4)
        else:
            self.y0 = arrset(self.y0,np.where(self.congrps==ig_['X2'])[0],float(4.4))
        if('X3' in ix_):
            self.x0[ix_['X3']] = float(12.02)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X3']),float(12.02)))
        if('X4' in ix_):
            self.x0[ix_['X4']] = float(11.82)
        else:
            self.y0 = arrset(self.y0,np.where(self.congrps==ig_['X4'])[0],float(11.82))
        if('X5' in ix_):
            self.x0[ix_['X5']] = float(0.702)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X5']),float(0.702)))
        if('X6' in ix_):
            self.x0[ix_['X6']] = float(0.852)
        else:
            self.y0 = arrset(self.y0,np.where(self.congrps==ig_['X6'])[0],float(0.852))
        pass
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eOE1', iet_)
        elftv = loaset(elftv,it,0,'X1')
        elftv = loaset(elftv,it,1,'X2')
        elftv = loaset(elftv,it,2,'X3')
        elftv = loaset(elftv,it,3,'X4')
        [it,iet_,_] = s2mpj_ii( 'eOE2', iet_)
        elftv = loaset(elftv,it,0,'X1')
        elftv = loaset(elftv,it,1,'X2')
        elftv = loaset(elftv,it,2,'X3')
        elftv = loaset(elftv,it,3,'X4')
        [it,iet_,_] = s2mpj_ii( 'eOE3', iet_)
        elftv = loaset(elftv,it,0,'X1')
        elftv = loaset(elftv,it,1,'X2')
        elftv = loaset(elftv,it,2,'X3')
        elftv = loaset(elftv,it,3,'X4')
        elftv = loaset(elftv,it,4,'X5')
        [it,iet_,_] = s2mpj_ii( 'eOE4', iet_)
        elftv = loaset(elftv,it,0,'X1')
        elftv = loaset(elftv,it,1,'X2')
        elftv = loaset(elftv,it,2,'X3')
        elftv = loaset(elftv,it,3,'X4')
        elftv = loaset(elftv,it,4,'X6')
        [it,iet_,_] = s2mpj_ii( 'eC1E1', iet_)
        elftv = loaset(elftv,it,0,'X1')
        elftv = loaset(elftv,it,1,'X2')
        elftv = loaset(elftv,it,2,'X3')
        elftv = loaset(elftv,it,3,'X4')
        elftv = loaset(elftv,it,4,'X5')
        elftv = loaset(elftv,it,5,'X6')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        ename = 'OE1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eOE1')
        ielftype = arrset(ielftype,ie,iet_["eOE1"])
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
        ename = 'OE2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eOE2')
        ielftype = arrset(ielftype,ie,iet_["eOE2"])
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
        ename = 'OE3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eOE3')
        ielftype = arrset(ielftype,ie,iet_["eOE3"])
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
        ename = 'OE4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eOE4')
        ielftype = arrset(ielftype,ie,iet_["eOE4"])
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
        vname = 'X6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X6')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'C1E1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eC1E1')
        ielftype = arrset(ielftype,ie,iet_["eC1E1"])
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
        vname = 'X6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X6')[0]
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
        self.grelt = loaset(self.grelt,ig,posel,ie_['OE1'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(2.04e-2))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['OE2'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.87e-2))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['OE3'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(6.07e-2))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['OE4'])
        self.grelw = loaset(self.grelw,ig,posel,float(4.37e-2))
        ig = ig_['C1']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['C1E1'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0e-3))
        ig = ig_['C2']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['OE3'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(6.2e-4))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['OE4'])
        self.grelw = loaset(self.grelw,ig,posel,float(5.8e-4))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               135.075961
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.cupper[np.arange(self.nle)] = np.zeros((self.nle,1))
        self.clower[np.arange(self.nle+self.neq,self.m)] = np.zeros((self.nge,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-COOR2-MY-6-2"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eOE1(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((3,4))
        IV_ = np.zeros(3)
        U_[0,0] = U_[0,0]+1
        U_[1,3] = U_[1,3]+1
        U_[2,0] = U_[2,0]+1
        U_[2,1] = U_[2,1]+1
        U_[2,2] = U_[2,2]+1
        IV_[0] = to_scalar(U_[0:1,:].dot(EV_))
        IV_[1] = to_scalar(U_[1:2,:].dot(EV_))
        IV_[2] = to_scalar(U_[2:3,:].dot(EV_))
        f_   = IV_[0]*IV_[1]*IV_[2]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = IV_[1]*IV_[2]
            g_[1] = IV_[0]*IV_[2]
            g_[2] = IV_[0]*IV_[1]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = IV_[2]
                H_[1,0] = H_[0,1]
                H_[0,2] = IV_[1]
                H_[2,0] = H_[0,2]
                H_[1,2] = IV_[0]
                H_[2,1] = H_[1,2]
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eOE2(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((3,4))
        IV_ = np.zeros(3)
        U_[0,1] = U_[0,1]+1
        U_[1,2] = U_[1,2]+1
        U_[2,0] = U_[2,0]+1
        U_[2,1] = U_[2,1]+1.570000e+00
        U_[2,3] = U_[2,3]+1
        IV_[0] = to_scalar(U_[0:1,:].dot(EV_))
        IV_[1] = to_scalar(U_[1:2,:].dot(EV_))
        IV_[2] = to_scalar(U_[2:3,:].dot(EV_))
        f_   = IV_[0]*IV_[1]*IV_[2]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = IV_[1]*IV_[2]
            g_[1] = IV_[0]*IV_[2]
            g_[2] = IV_[0]*IV_[1]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = IV_[2]
                H_[1,0] = H_[0,1]
                H_[0,2] = IV_[1]
                H_[2,0] = H_[0,2]
                H_[1,2] = IV_[0]
                H_[2,1] = H_[1,2]
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eOE3(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((4,5))
        IV_ = np.zeros(4)
        U_[0,0] = U_[0,0]+1
        U_[1,3] = U_[1,3]+1
        U_[2,4] = U_[2,4]+1
        U_[3,0] = U_[3,0]+1
        U_[3,1] = U_[3,1]+1
        U_[3,2] = U_[3,2]+1
        IV_[0] = to_scalar(U_[0:1,:].dot(EV_))
        IV_[1] = to_scalar(U_[1:2,:].dot(EV_))
        IV_[2] = to_scalar(U_[2:3,:].dot(EV_))
        IV_[3] = to_scalar(U_[3:4,:].dot(EV_))
        f_   = IV_[0]*IV_[1]*(IV_[2]**2)*IV_[3]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = IV_[1]*(IV_[2]**2)*IV_[3]
            g_[1] = IV_[0]*(IV_[2]**2)*IV_[3]
            g_[2] = IV_[0]*IV_[1]*2.0e+0*IV_[2]*IV_[3]
            g_[3] = IV_[0]*IV_[1]*(IV_[2]**2)
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,1] = (IV_[2]**2)*IV_[3]
                H_[1,0] = H_[0,1]
                H_[0,2] = IV_[1]*2.0e+0*IV_[2]*IV_[3]
                H_[2,0] = H_[0,2]
                H_[0,3] = IV_[1]*(IV_[2]**2)
                H_[3,0] = H_[0,3]
                H_[1,2] = IV_[0]*2.0e+0*IV_[2]*IV_[3]
                H_[2,1] = H_[1,2]
                H_[1,3] = IV_[0]*(IV_[2]**2)
                H_[3,1] = H_[1,3]
                H_[2,2] = IV_[0]*IV_[1]*2.0e+0*IV_[3]
                H_[2,3] = IV_[0]*IV_[1]*2.0e+0*IV_[2]
                H_[3,2] = H_[2,3]
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eOE4(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((4,5))
        IV_ = np.zeros(4)
        U_[0,1] = U_[0,1]+1
        U_[1,2] = U_[1,2]+1
        U_[2,4] = U_[2,4]+1
        U_[3,0] = U_[3,0]+1
        U_[3,1] = U_[3,1]+1.570000e+00
        U_[3,3] = U_[3,3]+1
        IV_[0] = to_scalar(U_[0:1,:].dot(EV_))
        IV_[1] = to_scalar(U_[1:2,:].dot(EV_))
        IV_[2] = to_scalar(U_[2:3,:].dot(EV_))
        IV_[3] = to_scalar(U_[3:4,:].dot(EV_))
        f_   = IV_[0]*IV_[1]*(IV_[2]**2)*IV_[3]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = IV_[1]*(IV_[2]**2)*IV_[3]
            g_[1] = IV_[0]*(IV_[2]**2)*IV_[3]
            g_[2] = IV_[0]*IV_[1]*2.0e+0*IV_[2]*IV_[3]
            g_[3] = IV_[0]*IV_[1]*(IV_[2]**2)
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,1] = (IV_[2]**2)*IV_[3]
                H_[1,0] = H_[0,1]
                H_[0,2] = IV_[1]*2.0e+0*IV_[2]*IV_[3]
                H_[2,0] = H_[0,2]
                H_[0,3] = IV_[1]*(IV_[2]**2)
                H_[3,0] = H_[0,3]
                H_[1,2] = IV_[0]*2.0e+0*IV_[2]*IV_[3]
                H_[2,1] = H_[1,2]
                H_[1,3] = IV_[0]*(IV_[2]**2)
                H_[3,1] = H_[1,3]
                H_[2,2] = IV_[0]*IV_[1]*2.0e+0*IV_[3]
                H_[2,3] = IV_[0]*IV_[1]*2.0e+0*IV_[2]
                H_[3,2] = H_[2,3]
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eC1E1(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]
            g_[1] = EV_[0,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]
            g_[2] = EV_[0,0]*EV_[1,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]
            g_[3] = EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[4,0]*EV_[5,0]
            g_[4] = EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[5,0]
            g_[5] = EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]
            if nargout>2:
                H_ = np.zeros((6,6))
                H_[0,1] = EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]
                H_[2,0] = H_[0,2]
                H_[0,3] = EV_[1,0]*EV_[2,0]*EV_[4,0]*EV_[5,0]
                H_[3,0] = H_[0,3]
                H_[0,4] = EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[5,0]
                H_[4,0] = H_[0,4]
                H_[0,5] = EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]
                H_[5,0] = H_[0,5]
                H_[1,2] = EV_[0,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]
                H_[2,1] = H_[1,2]
                H_[1,3] = EV_[0,0]*EV_[2,0]*EV_[4,0]*EV_[5,0]
                H_[3,1] = H_[1,3]
                H_[1,4] = EV_[0,0]*EV_[2,0]*EV_[3,0]*EV_[5,0]
                H_[4,1] = H_[1,4]
                H_[1,5] = EV_[0,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]
                H_[5,1] = H_[1,5]
                H_[2,3] = EV_[0,0]*EV_[1,0]*EV_[4,0]*EV_[5,0]
                H_[3,2] = H_[2,3]
                H_[2,4] = EV_[0,0]*EV_[1,0]*EV_[3,0]*EV_[5,0]
                H_[4,2] = H_[2,4]
                H_[2,5] = EV_[0,0]*EV_[1,0]*EV_[3,0]*EV_[4,0]
                H_[5,2] = H_[2,5]
                H_[3,4] = EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[5,0]
                H_[4,3] = H_[3,4]
                H_[3,5] = EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[4,0]
                H_[5,3] = H_[3,5]
                H_[4,5] = EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]
                H_[5,4] = H_[4,5]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

