from s2mpjlib import *
class  GROWTH(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : GROWTH
#    *********
#    GROWTH problem in 3 variables
# 
#    Fit the observed growth g(n) from Gaussian Elimination
#    with complete pivoting to a function of the form
#         U1 * n ** ( U2 + LOG(n) * U3 )
# 
#    SIF input: Nick Gould, Nov, 1991.
# 
#    classification = "NOR2-AN-3-12"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'GROWTH'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'GROWTH'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 3
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
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
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'G8')
        [ig,ig_,_] = s2mpj_ii('G9',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'G9')
        [ig,ig_,_] = s2mpj_ii('G10',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'G10')
        [ig,ig_,_] = s2mpj_ii('G11',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'G11')
        [ig,ig_,_] = s2mpj_ii('G12',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'G12')
        [ig,ig_,_] = s2mpj_ii('G13',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'G13')
        [ig,ig_,_] = s2mpj_ii('G14',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'G14')
        [ig,ig_,_] = s2mpj_ii('G15',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'G15')
        [ig,ig_,_] = s2mpj_ii('G16',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'G16')
        [ig,ig_,_] = s2mpj_ii('G18',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'G18')
        [ig,ig_,_] = s2mpj_ii('G20',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'G20')
        [ig,ig_,_] = s2mpj_ii('G25',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'G25')
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
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
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
        ig = ig_['G8']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['G8'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['G9']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['G9'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['G10']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['G10'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['G11']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['G11'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['G12']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['G12'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['G13']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['G13'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['G14']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['G14'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['G15']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['G15'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['G16']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['G16'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['G18']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['G18'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['G20']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['G20'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['G25']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['G25'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN                0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "NOR2-AN-3-12"
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

