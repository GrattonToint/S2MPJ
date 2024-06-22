from s2mpjlib import *
class  GROWTHLS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : GROWTHLS
#    *********
#    GROWTH problem in 3 variables
# 
#    Fit the observed growth g(n) from Gaussian Elimination
#    with complete pivoting to a function of the form
#         U1 * n ** ( U2 + LOG(n) * U3 )
# 
#    SIF input: Nick Gould, Nov, 1991, modified by Ph. Toint, March 1994.
# 
#    classification = "SUR2-AN-3-0"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'GROWTHLS'

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
        v_['N'] = 3
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        pb.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2mpj_ii('U1',ix_)
        pb.xnames=arrset(pb.xnames,iv,'U1')
        [iv,ix_,_] = s2mpj_ii('U2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'U2')
        [iv,ix_,_] = s2mpj_ii('U3',ix_)
        pb.xnames=arrset(pb.xnames,iv,'U3')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('G8',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('G9',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('G10',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('G11',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('G12',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('G13',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('G14',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('G15',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('G16',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('G18',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('G20',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('G25',ig_)
        gtype = arrset(gtype,ig,'<>')
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        pbm.gconst = arrset(pbm.gconst,ig_['G8'],float(8.0))
        pbm.gconst = arrset(pbm.gconst,ig_['G9'],float(8.4305))
        pbm.gconst = arrset(pbm.gconst,ig_['G10'],float(9.5294))
        pbm.gconst = arrset(pbm.gconst,ig_['G11'],float(10.4627))
        pbm.gconst = arrset(pbm.gconst,ig_['G12'],float(12.0))
        pbm.gconst = arrset(pbm.gconst,ig_['G13'],float(13.0205))
        pbm.gconst = arrset(pbm.gconst,ig_['G14'],float(14.5949))
        pbm.gconst = arrset(pbm.gconst,ig_['G15'],float(16.1078))
        pbm.gconst = arrset(pbm.gconst,ig_['G16'],float(18.0596))
        pbm.gconst = arrset(pbm.gconst,ig_['G18'],float(20.4569))
        pbm.gconst = arrset(pbm.gconst,ig_['G20'],float(24.25))
        pbm.gconst = arrset(pbm.gconst,ig_['G25'],float(32.9863))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower = np.zeros((pb.n,1))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.x0[ix_['U1']] = float(100.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eFIT', iet_)
        elftv = loaset(elftv,it,0,'U1')
        elftv = loaset(elftv,it,1,'U2')
        elftv = loaset(elftv,it,2,'U3')
        elftp = []
        elftp = loaset(elftp,it,0,'RN')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        ename = 'G8'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eFIT')
            ielftype = arrset( ielftype,ie,iet_['eFIT'])
        vname = 'U1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'U2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'U3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='RN')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(8.0))
        ename = 'G9'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eFIT')
            ielftype = arrset( ielftype,ie,iet_['eFIT'])
        vname = 'U1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'U2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'U3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='RN')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(9.0))
        ename = 'G10'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eFIT')
            ielftype = arrset( ielftype,ie,iet_['eFIT'])
        vname = 'U1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'U2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'U3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='RN')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(10.0))
        ename = 'G11'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eFIT')
            ielftype = arrset( ielftype,ie,iet_['eFIT'])
        vname = 'U1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'U2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'U3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='RN')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(11.0))
        ename = 'G12'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eFIT')
            ielftype = arrset( ielftype,ie,iet_['eFIT'])
        vname = 'U1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'U2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'U3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='RN')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(12.0))
        ename = 'G13'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eFIT')
            ielftype = arrset( ielftype,ie,iet_['eFIT'])
        vname = 'U1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'U2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'U3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='RN')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(13.0))
        ename = 'G14'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eFIT')
            ielftype = arrset( ielftype,ie,iet_['eFIT'])
        vname = 'U1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'U2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'U3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='RN')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(14.0))
        ename = 'G15'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eFIT')
            ielftype = arrset( ielftype,ie,iet_['eFIT'])
        vname = 'U1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'U2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'U3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='RN')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(15.0))
        ename = 'G16'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eFIT')
            ielftype = arrset( ielftype,ie,iet_['eFIT'])
        vname = 'U1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'U2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'U3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='RN')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(16.0))
        ename = 'G18'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eFIT')
            ielftype = arrset( ielftype,ie,iet_['eFIT'])
        vname = 'U1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'U2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'U3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='RN')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(18.0))
        ename = 'G20'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eFIT')
            ielftype = arrset( ielftype,ie,iet_['eFIT'])
        vname = 'U1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'U2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'U3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='RN')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(20.0))
        ename = 'G25'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eFIT')
            ielftype = arrset( ielftype,ie,iet_['eFIT'])
        vname = 'U1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'U2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'U3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='RN')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(25.0))
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
        ig = ig_['G8']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['G8'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['G9']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['G9'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['G10']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['G10'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['G11']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['G11'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['G12']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['G12'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['G13']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['G13'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['G14']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['G14'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['G15']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['G15'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['G16']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['G16'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['G18']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['G18'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['G20']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['G20'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['G25']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['G25'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN                0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "SUR2-AN-3-0"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eFIT(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        LOGRN = np.log(pbm.elpar[iel_][0])
        POWER = pbm.elpar[iel_][0]**(EV_[1]+LOGRN*EV_[2])
        f_   = EV_[0]*POWER
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = POWER
            g_[1] = EV_[0]*POWER*LOGRN
            g_[2] = EV_[0]*POWER*LOGRN**2
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = 0.0
                H_[0,1] = POWER*LOGRN
                H_[1,0] = H_[0,1]
                H_[0,2] = POWER*LOGRN**2
                H_[2,0] = H_[0,2]
                H_[1,1] = EV_[0]*POWER*LOGRN**2
                H_[1,2] = EV_[0]*POWER*LOGRN**3
                H_[2,1] = H_[1,2]
                H_[2,2] = EV_[0]*POWER*LOGRN**4
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

