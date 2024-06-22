from s2mpjlib import *
class  RES(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : RES
#    *********
# 
#    Dassault France ressort (spring) problem
# 
#    SIF input:  A. R. Conn, June 1993.
# 
#    classification = "NLR2-MN-20-14"
# 
# 
# 
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'RES'

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
        [iv,ix_,_] = s2mpj_ii('L0',ix_)
        pb.xnames=arrset(pb.xnames,iv,'L0')
        [iv,ix_,_] = s2mpj_ii('N',ix_)
        pb.xnames=arrset(pb.xnames,iv,'N')
        [iv,ix_,_] = s2mpj_ii('F',ix_)
        pb.xnames=arrset(pb.xnames,iv,'F')
        [iv,ix_,_] = s2mpj_ii('K',ix_)
        pb.xnames=arrset(pb.xnames,iv,'K')
        [iv,ix_,_] = s2mpj_ii('LB',ix_)
        pb.xnames=arrset(pb.xnames,iv,'LB')
        [iv,ix_,_] = s2mpj_ii('L',ix_)
        pb.xnames=arrset(pb.xnames,iv,'L')
        [iv,ix_,_] = s2mpj_ii('DE',ix_)
        pb.xnames=arrset(pb.xnames,iv,'DE')
        [iv,ix_,_] = s2mpj_ii('DI',ix_)
        pb.xnames=arrset(pb.xnames,iv,'DI')
        [iv,ix_,_] = s2mpj_ii('TO',ix_)
        pb.xnames=arrset(pb.xnames,iv,'TO')
        [iv,ix_,_] = s2mpj_ii('TOB',ix_)
        pb.xnames=arrset(pb.xnames,iv,'TOB')
        [iv,ix_,_] = s2mpj_ii('NU',ix_)
        pb.xnames=arrset(pb.xnames,iv,'NU')
        [iv,ix_,_] = s2mpj_ii('D',ix_)
        pb.xnames=arrset(pb.xnames,iv,'D')
        [iv,ix_,_] = s2mpj_ii('P',ix_)
        pb.xnames=arrset(pb.xnames,iv,'P')
        [iv,ix_,_] = s2mpj_ii('E',ix_)
        pb.xnames=arrset(pb.xnames,iv,'E')
        [iv,ix_,_] = s2mpj_ii('P0',ix_)
        pb.xnames=arrset(pb.xnames,iv,'P0')
        [iv,ix_,_] = s2mpj_ii('G',ix_)
        pb.xnames=arrset(pb.xnames,iv,'G')
        [iv,ix_,_] = s2mpj_ii('DM',ix_)
        pb.xnames=arrset(pb.xnames,iv,'DM')
        [iv,ix_,_] = s2mpj_ii('FR',ix_)
        pb.xnames=arrset(pb.xnames,iv,'FR')
        [iv,ix_,_] = s2mpj_ii('TOLIM',ix_)
        pb.xnames=arrset(pb.xnames,iv,'TOLIM')
        [iv,ix_,_] = s2mpj_ii('TOBLIM',ix_)
        pb.xnames=arrset(pb.xnames,iv,'TOBLIM')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('E1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E1')
        iv = ix_['F']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('E2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E2')
        iv = ix_['K']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('E3',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E3')
        iv = ix_['DE']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['D']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['DM']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('E4',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E4')
        iv = ix_['DI']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['D']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['DM']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('E5',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E5')
        iv = ix_['D']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['P']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['E']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('E6',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E6')
        iv = ix_['NU']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['N']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('E7',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E7')
        iv = ix_['D']
        pbm.A[ig,iv] = float(1.5)+pbm.A[ig,iv]
        iv = ix_['L0']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('E8',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E8')
        iv = ix_['L']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['LB']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['FR']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('E9',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E9')
        iv = ix_['LB']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('E10',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E10')
        iv = ix_['L']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['L0']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['F']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('E11',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E11')
        iv = ix_['TO']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('E12',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E12')
        iv = ix_['TOB']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('E13',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'E13')
        iv = ix_['TO']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['TOLIM']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('E14',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'E14')
        iv = ix_['TOB']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['TOBLIM']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        legrps = find(gtype,lambda x:x=='<=')
        eqgrps = find(gtype,lambda x:x=='==')
        gegrps = find(gtype,lambda x:x=='>=')
        pb.nle = len(legrps)
        pb.neq = len(eqgrps)
        pb.nge = len(gegrps)
        pb.m   = pb.nle+pb.neq+pb.nge
        pbm.congrps = find(gtype,lambda x:(x=='<=' or x=='==' or x=='>='))
        pb.cnames= cnames[pbm.congrps]
        pb.nob = ngrp-pb.m
        pbm.objgrps = find(gtype,lambda x:x=='<>')
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        pbm.gconst = arrset(pbm.gconst,ig_['E6'],float(-2.0))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),float('inf'))
        pb.xupper[ix_['L0']] = 100.0
        pb.xupper[ix_['N']] = 100.0
        pb.xupper[ix_['F']] = 30.0
        pb.xupper[ix_['K']] = 100.0
        pb.xupper[ix_['LB']] = 50.0
        pb.xupper[ix_['L']] = 50.0
        pb.xupper[ix_['DE']] = 30.0
        pb.xupper[ix_['DI']] = 30.0
        pb.xupper[ix_['TO']] = 800.0
        pb.xupper[ix_['TOB']] = 800.0
        pb.xupper[ix_['NU']] = 50.0
        pb.xlower[ix_['NU']] = 0.5
        pb.xupper[ix_['D']] = 10.0
        pb.xlower[ix_['D']] = 0.1
        pb.xupper[ix_['P']] = 20.0
        pb.xupper[ix_['E']] = 10.0
        pb.xupper[ix_['P0']] = 1000.0
        pb.xlower[ix_['P0']] = 1.0
        pb.xupper[ix_['G']] = 80000.0
        pb.xlower[ix_['G']] = 40000.0
        pb.xupper[ix_['DM']] = 30.0
        pb.xlower[ix_['DM']] = 0.1
        pb.xupper[ix_['FR']] = 50.0
        pb.xupper[ix_['TOLIM']] = 1000.0
        pb.xlower[ix_['TOLIM']] = 100.0
        pb.xupper[ix_['TOBLIM']] = 1000.0
        pb.xlower[ix_['TOBLIM']] = 100.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        pb.x0[ix_['L0']] = float(1.5000e-01)
        pb.x0[ix_['N']] = float(2.4079e+01)
        pb.x0[ix_['F']] = float(9.2459e-15)
        pb.x0[ix_['K']] = float(0.0000)
        pb.x0[ix_['LB']] = float(0.0000)
        pb.x0[ix_['L']] = float(1.5000e-01)
        pb.x0[ix_['DE']] = float(6.8120)
        pb.x0[ix_['DI']] = float(6.6120)
        pb.x0[ix_['TO']] = float(0.0000)
        pb.x0[ix_['TOB']] = float(0.0000)
        pb.x0[ix_['NU']] = float(2.2079e+01)
        pb.x0[ix_['D']] = float(1.0000e-01)
        pb.x0[ix_['P']] = float(6.5268e-01)
        pb.x0[ix_['E']] = float(5.5268e-01)
        pb.x0[ix_['P0']] = float(6.5887e+02)
        pb.x0[ix_['G']] = float(6.5887e+04)
        pb.x0[ix_['DM']] = float(6.7120)
        pb.x0[ix_['FR']] = float(1.5000e-01)
        pb.x0[ix_['TOLIM']] = float(1.0000e+02)
        pb.x0[ix_['TOBLIM']] = float(1.0000e+02)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        [it,iet_,_] = s2mpj_ii( 'en311d14', iet_)
        elftv = loaset(elftv,it,0,'V')
        elftv = loaset(elftv,it,1,'W')
        elftv = loaset(elftv,it,2,'X')
        elftv = loaset(elftv,it,3,'Y')
        elftv = loaset(elftv,it,4,'Z')
        [it,iet_,_] = s2mpj_ii( 'en14d31', iet_)
        elftv = loaset(elftv,it,0,'W')
        elftv = loaset(elftv,it,1,'X')
        elftv = loaset(elftv,it,2,'Y')
        elftv = loaset(elftv,it,3,'Z')
        [it,iet_,_] = s2mpj_ii( 'en11d3', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        elftv = loaset(elftv,it,2,'Z')
        [it,iet_,_] = s2mpj_ii( 'en111d2', iet_)
        elftv = loaset(elftv,it,0,'W')
        elftv = loaset(elftv,it,1,'X')
        elftv = loaset(elftv,it,2,'Y')
        elftv = loaset(elftv,it,3,'Z')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        ename = 'EL1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en311d14')
        ielftype = arrset(ielftype, ie, iet_["en311d14"])
        vname = 'DM'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NU'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='W')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'P0'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'G'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'D'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EL2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en14d31')
        ielftype = arrset(ielftype, ie, iet_["en14d31"])
        vname = 'G'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='W')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'D'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'DM'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NU'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EL3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
        ielftype = arrset(ielftype, ie, iet_["en2PR"])
        vname = 'NU'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'P'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EL4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en11d3')
        ielftype = arrset(ielftype, ie, iet_["en11d3"])
        vname = 'P0'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'DM'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'D'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EL5'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en111d2')
        ielftype = arrset(ielftype, ie, iet_["en111d2"])
        vname = 'G'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='W')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'D'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'E'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'DM'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.cupper[np.arange(pb.nle)] = np.zeros((pb.nle,1))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.lincons   = np.arange(len(pbm.congrps))
        pb.pbclass = "NLR2-MN-20-14"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def e_globs(pbm):

        import numpy as np
        pbm.efpar = np.array([]);
        pbm.efpar = arrset( pbm.efpar,0,3.1415926535)
        return pbm

    @staticmethod
    def en2PR(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]
            g_[1] = EV_[0]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = 1.0
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def en311d14(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        V3WX = EV_[0]**3*EV_[1]*EV_[2]
        YZ4 = EV_[3]*EV_[4]**4
        f_   = V3WX/YZ4
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = (3.0*EV_[0]**2*EV_[1]*EV_[2])/YZ4
            g_[1] = (EV_[0]**3*EV_[2])/YZ4
            g_[2] = (EV_[0]**3*EV_[1])/YZ4
            g_[3] = -V3WX/(EV_[3]*YZ4)
            g_[4] = -(4.0*V3WX)/(EV_[4]*YZ4)
            if nargout>2:
                H_ = np.zeros((5,5))
                H_[0,0] = (6.0*EV_[0]*EV_[1]*EV_[2])/YZ4
                H_[0,1] = (3.0*EV_[0]**2*EV_[2])/YZ4
                H_[1,0] = H_[0,1]
                H_[0,2] = (3.0*EV_[0]**2*EV_[1])/YZ4
                H_[2,0] = H_[0,2]
                H_[0,3] = -(3.0*EV_[0]**2*EV_[1])/(YZ4*EV_[3])
                H_[3,0] = H_[0,3]
                H_[0,4] = -(12.0*EV_[0]**2*EV_[1]*EV_[2])/(YZ4*EV_[4])
                H_[4,0] = H_[0,4]
                H_[1,2] = EV_[0]**3/YZ4
                H_[2,1] = H_[1,2]
                H_[1,3] = -(EV_[0]**3*EV_[2])/(YZ4*EV_[3])
                H_[3,1] = H_[1,3]
                H_[1,4] = -(4.0*EV_[0]**3*EV_[2])/(YZ4*EV_[4])
                H_[4,1] = H_[1,4]
                H_[2,3] = -(EV_[0]**3*EV_[1])/(YZ4*EV_[3])
                H_[3,2] = H_[2,3]
                H_[2,4] = -(4.0*EV_[0]**3*EV_[1])/(YZ4*EV_[4])
                H_[4,2] = H_[2,4]
                H_[3,3] = -(2.0*V3WX)/(EV_[3]**2*YZ4)
                H_[3,4] = (4.0*V3WX)/(EV_[3]*EV_[4]*YZ4)
                H_[4,3] = H_[3,4]
                H_[4,4] = (20.0*V3WX)/(EV_[4]**2*YZ4)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def en14d31(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        WX4 = EV_[0]*EV_[1]**4
        Y3Z = EV_[2]**3*EV_[3]
        f_   = WX4/Y3Z
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]**4/Y3Z
            g_[1] = (4.0*EV_[0]*EV_[1]**3)/Y3Z
            g_[2] = -(3.0*WX4)/(EV_[2]*Y3Z)
            g_[3] = -WX4/(EV_[3]*Y3Z)
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,1] = (4.0*EV_[1]**3)/Y3Z
                H_[1,0] = H_[0,1]
                H_[0,2] = -(3.0*EV_[1]**4)/(EV_[2]*Y3Z)
                H_[2,0] = H_[0,2]
                H_[0,3] = -EV_[1]**4/(EV_[3]*Y3Z)
                H_[3,0] = H_[0,3]
                H_[1,1] = (12.0*EV_[0]*EV_[1]**2)/Y3Z
                H_[1,2] = -(12.0*EV_[0]*EV_[1]**3)/(EV_[2]*Y3Z)
                H_[2,1] = H_[1,2]
                H_[1,3] = -(4.0*EV_[0]*EV_[1]**3)/(EV_[3]*Y3Z)
                H_[3,1] = H_[1,3]
                H_[2,2] = (12.0*WX4)/(EV_[2]**2*Y3Z)
                H_[2,3] = (3.0*WX4)/(EV_[3]*EV_[2]*Y3Z)
                H_[3,2] = H_[2,3]
                H_[3,3] = (2.0*WX4)/(EV_[3]**2*Y3Z)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def en11d3(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = (EV_[0]*EV_[1])/(pbm.efpar[0]*EV_[2]**3)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]/(pbm.efpar[0]*EV_[2]**3)
            g_[1] = EV_[0]/(pbm.efpar[0]*EV_[2]**3)
            g_[2] = -(3.0*EV_[0]*EV_[1])/(pbm.efpar[0]*EV_[2]**4)
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = 1.0/(pbm.efpar[0]*EV_[2]**3)
                H_[1,0] = H_[0,1]
                H_[0,2] = -(3.0*EV_[1])/(pbm.efpar[0]*EV_[2]**4)
                H_[2,0] = H_[0,2]
                H_[1,2] = -(3.0*EV_[0])/(pbm.efpar[0]*EV_[2]**4)
                H_[2,1] = H_[1,2]
                H_[2,2] = (12.0*EV_[0]*EV_[1])/(pbm.efpar[0]*EV_[2]**5)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def en111d2(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = (EV_[0]*EV_[1]*EV_[2])/(pbm.efpar[0]*EV_[3]**2)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = (EV_[1]*EV_[2])/(pbm.efpar[0]*EV_[3]**2)
            g_[1] = (EV_[0]*EV_[2])/(pbm.efpar[0]*EV_[3]**2)
            g_[2] = (EV_[0]*EV_[1])/(pbm.efpar[0]*EV_[3]**2)
            g_[3] = -(2.0*EV_[0]*EV_[1]*EV_[2])/(pbm.efpar[0]*EV_[3]**3)
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,1] = EV_[2]/(pbm.efpar[0]*EV_[3]**2)
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1]/(pbm.efpar[0]*EV_[3]**2)
                H_[2,0] = H_[0,2]
                H_[0,3] = -(2.0*EV_[1]*EV_[2])/(pbm.efpar[0]*EV_[3]**3)
                H_[3,0] = H_[0,3]
                H_[1,2] = EV_[0]/(pbm.efpar[0]*EV_[3]**2)
                H_[2,1] = H_[1,2]
                H_[1,3] = -(2.0*EV_[0]*EV_[2])/(pbm.efpar[0]*EV_[3]**3)
                H_[3,1] = H_[1,3]
                H_[2,3] = -(2.0*EV_[0]*EV_[1])/(pbm.efpar[0]*EV_[3]**3)
                H_[3,2] = H_[2,3]
                H_[3,3] = (6.0*EV_[0]*EV_[1]*EV_[2])/(pbm.efpar[0]*EV_[3]**4)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

