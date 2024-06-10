from s2mpjlib import *
class  SYNTHES2(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : SYNTHES2
#    *********
# 
#    Source: Test problem 2 (Synthesis of processing system) in 
#    M. Duran & I.E. Grossmann,
#    "An outer approximation algorithm for a class of mixed integer nonlinear
#     programs", Mathematical Programming 36, pp. 307-339, 1986.
# 
#    SIF input: S. Leyffer, October 1997
# 
#    classification = "OOR2-AN-11-14"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'SYNTHES2'

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
        v_['1'] = 1
        v_['5'] = 5
        v_['6'] = 6
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['6'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
        for I in range(int(v_['1']),int(v_['5'])+1):
            [iv,ix_,_] = s2mpj_ii('Y'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'Y'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['Y1']
        pbm.A[ig,iv] = float(5.0)+pbm.A[ig,iv]
        iv = ix_['Y2']
        pbm.A[ig,iv] = float(8.0)+pbm.A[ig,iv]
        iv = ix_['Y3']
        pbm.A[ig,iv] = float(6.0)+pbm.A[ig,iv]
        iv = ix_['Y4']
        pbm.A[ig,iv] = float(10.0)+pbm.A[ig,iv]
        iv = ix_['Y5']
        pbm.A[ig,iv] = float(6.0)+pbm.A[ig,iv]
        iv = ix_['X1']
        pbm.A[ig,iv] = float(-10.0)+pbm.A[ig,iv]
        iv = ix_['X2']
        pbm.A[ig,iv] = float(-15.0)+pbm.A[ig,iv]
        iv = ix_['X3']
        pbm.A[ig,iv] = float(-15.0)+pbm.A[ig,iv]
        iv = ix_['X4']
        pbm.A[ig,iv] = float(15.0)+pbm.A[ig,iv]
        iv = ix_['X5']
        pbm.A[ig,iv] = float(5.0)+pbm.A[ig,iv]
        iv = ix_['X6']
        pbm.A[ig,iv] = float(-20.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('N1',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'N1')
        [ig,ig_,_] = s2mpj_ii('N2',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'N2')
        iv = ix_['Y1']
        pbm.A[ig,iv] = float(-10.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('N3',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'N3')
        iv = ix_['Y2']
        pbm.A[ig,iv] = float(-10.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('L1',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'L1')
        iv = ix_['X3']
        pbm.A[ig,iv] = float(1.25)+pbm.A[ig,iv]
        iv = ix_['Y3']
        pbm.A[ig,iv] = float(-10.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('L2',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'L2')
        iv = ix_['X4']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X5']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['Y4']
        pbm.A[ig,iv] = float(-10.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('L3',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'L3')
        iv = ix_['X3']
        pbm.A[ig,iv] = float(-2.0)+pbm.A[ig,iv]
        iv = ix_['X6']
        pbm.A[ig,iv] = float(2.0)+pbm.A[ig,iv]
        iv = ix_['Y5']
        pbm.A[ig,iv] = float(-10.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('L4',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'L4')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X2']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X3']
        pbm.A[ig,iv] = float(-2.0)+pbm.A[ig,iv]
        iv = ix_['X4']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X6']
        pbm.A[ig,iv] = float(2.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('L5',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'L5')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X2']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X3']
        pbm.A[ig,iv] = float(-0.75)+pbm.A[ig,iv]
        iv = ix_['X4']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X6']
        pbm.A[ig,iv] = float(2.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('L6',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'L6')
        iv = ix_['X3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X6']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('L7',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'L7')
        iv = ix_['X3']
        pbm.A[ig,iv] = float(2.0)+pbm.A[ig,iv]
        iv = ix_['X4']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X6']
        pbm.A[ig,iv] = float(-2.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('L8',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'L8')
        iv = ix_['X4']
        pbm.A[ig,iv] = float(-0.5)+pbm.A[ig,iv]
        iv = ix_['X5']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('L9',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'L9')
        iv = ix_['X4']
        pbm.A[ig,iv] = float(0.2)+pbm.A[ig,iv]
        iv = ix_['X5']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('L10',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'L10')
        iv = ix_['Y1']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['Y2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('L11',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'L11')
        iv = ix_['Y4']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['Y5']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
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
        pbm.gconst = arrset(pbm.gconst,ig_['OBJ'],float(-140.0))
        pbm.gconst = arrset(pbm.gconst,ig_['N2'],float(1.0))
        pbm.gconst = arrset(pbm.gconst,ig_['N3'],float(1.0))
        pbm.gconst = arrset(pbm.gconst,ig_['L10'],float(1.0))
        pbm.gconst = arrset(pbm.gconst,ig_['L11'],float(1.0))
        #%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange = np.full((ngrp,1),None)
        grange[legrps] = np.full((pb.nle,1),float('inf'))
        grange = arrset(grange,ig_['L11'],float(-1.0))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),float('inf'))
        pb.xupper[ix_['X1']] = 2.0
        pb.xupper[ix_['X2']] = 2.0
        pb.xupper[ix_['X3']] = 2.0
        pb.xupper[ix_['X6']] = 3.0
        for I in range(int(v_['1']),int(v_['5'])+1):
            pb.xupper[ix_['Y'+str(I)]] = 1.0
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eLOGSUM', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        [it,iet_,_] = s2mpj_ii( 'eEXPA', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftp = []
        elftp = loaset(elftp,it,0,'A')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        ename = 'LOGX4X5'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eLOGSUM')
        ielftype = arrset(ielftype, ie, iet_["eLOGSUM"])
        pb.x0 = np.zeros((pb.n,1))
        vname = 'X4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EXPX1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eEXPA')
        ielftype = arrset(ielftype, ie, iet_["eEXPA"])
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='A')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        ename = 'EXPX2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eEXPA')
        ielftype = arrset(ielftype, ie, iet_["eEXPA"])
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='A')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.2))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        ig = ig_['OBJ']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EXPX1'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EXPX2'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['LOGX4X5'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-60.0))
        ig = ig_['N1']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['LOGX4X5'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        ig = ig_['N2']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EXPX1'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
        ig = ig_['N3']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EXPX2'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle)] = grange[legrps]
        pb.cupper[np.arange(pb.nle)] = np.zeros((pb.nle,1))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "OOR2-AN-11-14"
        pb.x0          = np.zeros((pb.n,1))
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eLOGSUM(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = np.log(EV_[0]+EV_[1]+1.0)
        if not isinstance( f_, float ):
            f_   = f_.item();
        DX = 1.0/(EV_[0]+EV_[1]+1.0)
        DXDX = -1.0/(EV_[0]+EV_[1]+1.0)**2
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = DX
            g_[1] = DX
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = DXDX
                H_[0,1] = DXDX
                H_[1,0] = H_[0,1]
                H_[1,1] = DXDX
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eEXPA(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        EXPXA = np.exp(EV_[0]/pbm.elpar[iel_][0])
        f_   = EXPXA
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EXPXA/pbm.elpar[iel_][0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = EXPXA/pbm.elpar[iel_][0]/pbm.elpar[iel_][0]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

