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
#    classification = "C-CSUR2-AN-3-0"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'GROWTHLS'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 3
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [iv,ix_,_] = s2mpj_ii('U1',ix_)
        self.xnames=arrset(self.xnames,iv,'U1')
        [iv,ix_,_] = s2mpj_ii('U2',ix_)
        self.xnames=arrset(self.xnames,iv,'U2')
        [iv,ix_,_] = s2mpj_ii('U3',ix_)
        self.xnames=arrset(self.xnames,iv,'U3')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
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
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        selfnob      = ngrp
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        self.gconst = arrset(self.gconst,ig_['G8'],float(8.0))
        self.gconst = arrset(self.gconst,ig_['G9'],float(8.4305))
        self.gconst = arrset(self.gconst,ig_['G10'],float(9.5294))
        self.gconst = arrset(self.gconst,ig_['G11'],float(10.4627))
        self.gconst = arrset(self.gconst,ig_['G12'],float(12.0))
        self.gconst = arrset(self.gconst,ig_['G13'],float(13.0205))
        self.gconst = arrset(self.gconst,ig_['G14'],float(14.5949))
        self.gconst = arrset(self.gconst,ig_['G15'],float(16.1078))
        self.gconst = arrset(self.gconst,ig_['G16'],float(18.0596))
        self.gconst = arrset(self.gconst,ig_['G18'],float(20.4569))
        self.gconst = arrset(self.gconst,ig_['G20'],float(24.25))
        self.gconst = arrset(self.gconst,ig_['G25'],float(32.9863))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.x0[ix_['U1']] = float(100.0)
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
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        ename = 'G8'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eFIT')
            ielftype = arrset(ielftype,ie,iet_['eFIT'])
        vname = 'U1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='RN')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(8.0))
        ename = 'G9'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eFIT')
            ielftype = arrset(ielftype,ie,iet_['eFIT'])
        vname = 'U1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='RN')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(9.0))
        ename = 'G10'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eFIT')
            ielftype = arrset(ielftype,ie,iet_['eFIT'])
        vname = 'U1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='RN')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(10.0))
        ename = 'G11'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eFIT')
            ielftype = arrset(ielftype,ie,iet_['eFIT'])
        vname = 'U1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='RN')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(11.0))
        ename = 'G12'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eFIT')
            ielftype = arrset(ielftype,ie,iet_['eFIT'])
        vname = 'U1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='RN')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(12.0))
        ename = 'G13'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eFIT')
            ielftype = arrset(ielftype,ie,iet_['eFIT'])
        vname = 'U1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='RN')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(13.0))
        ename = 'G14'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eFIT')
            ielftype = arrset(ielftype,ie,iet_['eFIT'])
        vname = 'U1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='RN')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(14.0))
        ename = 'G15'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eFIT')
            ielftype = arrset(ielftype,ie,iet_['eFIT'])
        vname = 'U1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='RN')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(15.0))
        ename = 'G16'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eFIT')
            ielftype = arrset(ielftype,ie,iet_['eFIT'])
        vname = 'U1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='RN')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(16.0))
        ename = 'G18'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eFIT')
            ielftype = arrset(ielftype,ie,iet_['eFIT'])
        vname = 'U1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='RN')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(18.0))
        ename = 'G20'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eFIT')
            ielftype = arrset(ielftype,ie,iet_['eFIT'])
        vname = 'U1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='RN')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(20.0))
        ename = 'G25'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eFIT')
            ielftype = arrset(ielftype,ie,iet_['eFIT'])
        vname = 'U1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='RN')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(25.0))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gL2',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for ig in range(0,ngrp):
            self.grftype = arrset(self.grftype,ig,'gL2')
        ig = ig_['G8']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['G8'])
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['G9']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['G9'])
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['G10']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['G10'])
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['G11']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['G11'])
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['G12']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['G12'])
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['G13']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['G13'])
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['G14']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['G14'])
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['G15']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['G15'])
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['G16']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['G16'])
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['G18']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['G18'])
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['G20']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['G20'])
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['G25']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['G25'])
        self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 0.0
#    Solution
# LO SOLTN                0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass   = "C-CSUR2-AN-3-0"
        self.objderlvl = 2

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eFIT(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        LOGRN = np.log(self.elpar[iel_][0])
        POWER = self.elpar[iel_][0]**(EV_[1,0]+LOGRN*EV_[2,0])
        f_   = EV_[0,0]*POWER
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = POWER
            g_[1] = EV_[0,0]*POWER*LOGRN
            g_[2] = EV_[0,0]*POWER*LOGRN**2
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = 0.0
                H_[0,1] = POWER*LOGRN
                H_[1,0] = H_[0,1]
                H_[0,2] = POWER*LOGRN**2
                H_[2,0] = H_[0,2]
                H_[1,1] = EV_[0,0]*POWER*LOGRN**2
                H_[1,2] = EV_[0,0]*POWER*LOGRN**3
                H_[2,1] = H_[1,2]
                H_[2,2] = EV_[0,0]*POWER*LOGRN**4
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gL2(self,nargout,*args):

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

