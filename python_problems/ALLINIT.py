from s2mpjlib import *
class  ALLINIT(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : ALLINIT
#    *********
# 
#    A problem with "all in it". Intended to verify that changes
#    to LANCELOT are safe.
# 
#    Source:
#    N. Gould, private communication.
# 
#    SIF input: Nick Gould, June 1990.
# 
#    classification = "OBR2-AY-4-0"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'ALLINIT'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2mpj_ii('X1',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X1')
        [iv,ix_,_] = s2mpj_ii('X2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X2')
        [iv,ix_,_] = s2mpj_ii('X3',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X3')
        [iv,ix_,_] = s2mpj_ii('X4',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X4')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('FT1',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('FT2',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('FT3',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('FT4',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('FT5',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X4']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('FT6',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('FNT1',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('FNT2',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X4']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('FNT3',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('FNT4',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('FNT5',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('FNT6',ig_)
        gtype = arrset(gtype,ig,'<>')
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        pbm.gconst = arrset(pbm.gconst,ig_['FT2'],float(1.0))
        pbm.gconst = arrset(pbm.gconst,ig_['FT5'],float(3.0))
        pbm.gconst = arrset(pbm.gconst,ig_['FNT2'],float(1.0))
        pbm.gconst = arrset(pbm.gconst,ig_['FNT5'],float(4.0))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),float('inf'))
        pb.xlower[ix_['X1']] = -float('Inf')
        pb.xupper[ix_['X1']] = +float('Inf')
        pb.xlower[ix_['X2']] = 1.0
        pb.xlower[ix_['X3']] = -1.0e+10
        pb.xupper[ix_['X3']] = 1.0
        pb.xlower[ix_['X4']] = 2.0
        pb.xupper[ix_['X4']] = 2.0
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQR', iet_)
        elftv = loaset(elftv,it,0,'X')
        [it,iet_,_] = s2mpj_ii( 'eSQR2', iet_)
        elftv = loaset(elftv,it,0,'Y')
        elftv = loaset(elftv,it,1,'Z')
        [it,iet_,_] = s2mpj_ii( 'eSINSQR', iet_)
        elftv = loaset(elftv,it,0,'X')
        [it,iet_,_] = s2mpj_ii( 'ePRODSQR', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        ename = 'FT3E1'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eSQR')
            ielftype = arrset( ielftype,ie,iet_['eSQR'])
        pb.x0 = np.zeros((pb.n,1))
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FT4E1'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eSQR')
            ielftype = arrset( ielftype,ie,iet_['eSQR'])
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FT4E2'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSQR2')
        ielftype = arrset(ielftype, ie, iet_["eSQR2"])
        vname = 'X3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FT56E1'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSINSQR')
        ielftype = arrset(ielftype, ie, iet_["eSINSQR"])
        vname = 'X3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FT5E2'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePRODSQR')
        ielftype = arrset(ielftype, ie, iet_["ePRODSQR"])
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FNT3E1'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eSQR')
            ielftype = arrset( ielftype,ie,iet_['eSQR'])
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FNT4E1'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eSQR')
            ielftype = arrset( ielftype,ie,iet_['eSQR'])
        vname = 'X3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FNT4E2'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSQR2')
        ielftype = arrset(ielftype, ie, iet_["eSQR2"])
        vname = 'X4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FNT56E1'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSINSQR')
        ielftype = arrset(ielftype, ie, iet_["eSINSQR"])
        vname = 'X4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'FNT5E2'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePRODSQR')
        ielftype = arrset(ielftype, ie, iet_["ePRODSQR"])
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gTRIVIAL',igt_)
        [it,igt_,_] = s2mpj_ii('gL2',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for ig in range(0,ngrp):
            pbm.grftype = arrset(pbm.grftype,ig,'gL2')
        ig = ig_['FT1']
        pbm.grftype = arrset(pbm.grftype,ig,'gTRIVIAL')
        ig = ig_['FT2']
        pbm.grftype = arrset(pbm.grftype,ig,'gTRIVIAL')
        ig = ig_['FT3']
        pbm.grftype = arrset(pbm.grftype,ig,'gTRIVIAL')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FT3E1'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['FT4']
        pbm.grftype = arrset(pbm.grftype,ig,'gTRIVIAL')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FT4E1'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FT4E2'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_['FT5']
        pbm.grftype = arrset(pbm.grftype,ig,'gTRIVIAL')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FT56E1'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FT5E2'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_['FT6']
        pbm.grftype = arrset(pbm.grftype,ig,'gTRIVIAL')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FT56E1'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['FNT3']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FNT3E1'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['FNT4']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FNT4E1'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FNT4E2'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_['FNT5']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FNT56E1'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FNT5E2'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_['FNT6']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['FNT56E1'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "OBR2-AY-4-0"
        pb.x0          = np.zeros((pb.n,1))
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSQR(pbm,nargout,*args):

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
    def eSQR2(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((1,2))
        IV_ = np.zeros(1)
        U_[0,0] = U_[0,0]+1
        U_[0,1] = U_[0,1]+1
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

    @staticmethod
    def eSINSQR(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        SINX = np.sin(EV_[0])
        COSX = np.cos(EV_[0])
        f_   = SINX*SINX
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*SINX*COSX
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0*(COSX*COSX-SINX*SINX)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ePRODSQR(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        XX = EV_[0]*EV_[0]
        YY = EV_[1]*EV_[1]
        f_   = XX*YY
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*EV_[0]*YY
            g_[1] = 2.0*XX*EV_[1]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 2.0*YY
                H_[0,1] = 4.0*EV_[0]*EV_[1]
                H_[1,0] = H_[0,1]
                H_[1,1] = 2.0*XX
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gTRIVIAL(pbm,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= GVAR_
        if nargout>1:
            g_ = 1.0
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 0.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def gL2(pbm,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= GVAR_*GVAR_
        if nargout>1:
            g_ = GVAR_+GVAR_
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

