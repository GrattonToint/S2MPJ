from s2mpjlib import *
class  ENGVAL2(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : ENGVAL2
#    *********
#    The ENGVAL2 problem.
# 
#    Source: problem 15 in
#    A.R. Buckley,
#    "Test functions for unconstrained minimization",
#    TR 1989CS-3, Mathematics, statistics and computing centre,
#    Dalhousie University, Halifax (CDN), 1989.
# 
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "SUR2-AN-3-0"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'ENGVAL2'

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
        pb.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2mpj_ii('X1',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X1')
        [iv,ix_,_] = s2mpj_ii('X2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X2')
        [iv,ix_,_] = s2mpj_ii('X3',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X3')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('G1',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('G2',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('G3',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('G4',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X3']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('G5',ig_)
        gtype = arrset(gtype,ig,'<>')
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        pbm.gconst = arrset(pbm.gconst,ig_['G1'],float(1.0))
        pbm.gconst = arrset(pbm.gconst,ig_['G2'],float(1.0))
        pbm.gconst = arrset(pbm.gconst,ig_['G3'],float(1.0))
        pbm.gconst = arrset(pbm.gconst,ig_['G4'],float(-1.0))
        pbm.gconst = arrset(pbm.gconst,ig_['G5'],float(36.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower = np.zeros((pb.n,1))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.x0[ix_['X1']] = float(1.0)
        pb.x0[ix_['X2']] = float(2.0)
        pb.x0[ix_['X3']] = float(0.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'X')
        [it,iet_,_] = s2mpj_ii( 'eSSQ', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftp = []
        elftp = loaset(elftp,it,0,'S')
        [it,iet_,_] = s2mpj_ii( 'eELT', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        ename = 'E1'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
            ielftype = arrset( ielftype,ie,iet_['eSQ'])
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E2'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
            ielftype = arrset( ielftype,ie,iet_['eSQ'])
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E3'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
            ielftype = arrset( ielftype,ie,iet_['eSQ'])
        vname = 'X3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E4'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
            ielftype = arrset( ielftype,ie,iet_['eSQ'])
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E5'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eELT')
        ielftype = arrset(ielftype, ie, iet_["eELT"])
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E6'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSSQ')
        ielftype = arrset(ielftype, ie, iet_["eSSQ"])
        vname = 'X3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='S')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(2.0))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
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
        ig = ig_['G1']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E1'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E2'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E3'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['G2']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E1'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E2'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E6'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['G5']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E4'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(3.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E5'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN               0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "SUR2-AN-3-0"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSQ(pbm,nargout,*args):

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
    def eSSQ(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = (EV_[0]-pbm.elpar[iel_][0])**2
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*(EV_[0]-pbm.elpar[iel_][0])
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
    def eELT(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        TWOW = 5.0*EV_[1]-EV_[0]+1.0
        f_   = EV_[0]**3+TWOW*TWOW
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 3.0*EV_[0]**2-2.0*TWOW
            g_[1] = 10.0*TWOW
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 6.0*EV_[0]+2.0
                H_[0,1] = -10.0
                H_[1,0] = H_[0,1]
                H_[1,1] = 50.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

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

