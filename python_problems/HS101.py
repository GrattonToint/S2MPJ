from s2mpjlib import *
class  HS101(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    Source: problem 101 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: N. Gould, December 1989.
# 
#    classification = "C-COOR2-AN-7-5"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS101'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['1'] = 1
        v_['M'] = 5
        v_['N'] = 7
        v_['A101'] = -0.25
        v_['A102'] = 0.125
        v_['A103'] = 0.5
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
        for I in range(int(v_['1']),int(v_['M'])+1):
            [ig,ig_,_] = s2mpj_ii('CONSTR'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'CONSTR'+str(I))
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
        self.gconst = arrset(self.gconst,ig_['CONSTR1'],float(1.0))
        self.gconst = arrset(self.gconst,ig_['CONSTR2'],float(1.0))
        self.gconst = arrset(self.gconst,ig_['CONSTR3'],float(1.0))
        self.gconst = arrset(self.gconst,ig_['CONSTR4'],float(1.0))
        self.gconst = arrset(self.gconst,ig_['CONSTR5'],float(3000.0))
        #%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange = np.full((ngrp,1),None)
        grange[legrps] = -np.full((self.nle,1),float('inf'))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),0.1)
        self.xupper = np.full((self.n,1),10.0)
        self.xlower[ix_['X7']] = 0.01
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(6.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'en3PR', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftp = []
        elftp = loaset(elftp,it,0,'P1')
        elftp = loaset(elftp,it,1,'P2')
        elftp = loaset(elftp,it,2,'P3')
        [it,iet_,_] = s2mpj_ii( 'en4PR', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftv = loaset(elftv,it,3,'V4')
        elftp = loaset(elftp,it,0,'P1')
        elftp = loaset(elftp,it,1,'P2')
        elftp = loaset(elftp,it,2,'P3')
        elftp = loaset(elftp,it,3,'P4')
        [it,iet_,_] = s2mpj_ii( 'en5PR', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftv = loaset(elftv,it,3,'V4')
        elftv = loaset(elftv,it,4,'V5')
        elftp = loaset(elftp,it,0,'P1')
        elftp = loaset(elftp,it,1,'P2')
        elftp = loaset(elftp,it,2,'P3')
        elftp = loaset(elftp,it,3,'P4')
        elftp = loaset(elftp,it,4,'P5')
        [it,iet_,_] = s2mpj_ii( 'en6PR', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftv = loaset(elftv,it,3,'V4')
        elftv = loaset(elftv,it,4,'V5')
        elftv = loaset(elftv,it,5,'V6')
        elftp = loaset(elftp,it,0,'P1')
        elftp = loaset(elftp,it,1,'P6')
        elftp = loaset(elftp,it,2,'P2')
        elftp = loaset(elftp,it,3,'P3')
        elftp = loaset(elftp,it,4,'P4')
        elftp = loaset(elftp,it,5,'P5')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        ename = 'E1C1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en4PR')
        ielftype = arrset(ielftype,ie,iet_["en4PR"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V4')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='P1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(0.5))
        posep = np.where(elftp[ielftype[ie]]=='P2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.0))
        posep = np.where(elftp[ielftype[ie]]=='P3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-2.0))
        posep = np.where(elftp[ielftype[ie]]=='P4')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
        ename = 'E2C1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en5PR')
        ielftype = arrset(ielftype,ie,iet_["en5PR"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V4')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V5')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='P1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(3.0))
        posep = np.where(elftp[ielftype[ie]]=='P2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
        posep = np.where(elftp[ielftype[ie]]=='P3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-2.0))
        posep = np.where(elftp[ielftype[ie]]=='P4')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
        posep = np.where(elftp[ielftype[ie]]=='P5')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(0.5))
        ename = 'E3C1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en5PR')
        ielftype = arrset(ielftype,ie,iet_["en5PR"])
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V4')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V5')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='P1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.0))
        posep = np.where(elftp[ielftype[ie]]=='P2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
        posep = np.where(elftp[ielftype[ie]]=='P3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-0.5))
        posep = np.where(elftp[ielftype[ie]]=='P4')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(0.66666666))
        posep = np.where(elftp[ielftype[ie]]=='P5')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(0.25))
        ename = 'E1C2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en5PR')
        ielftype = arrset(ielftype,ie,iet_["en5PR"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V4')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V5')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='P1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-0.5))
        posep = np.where(elftp[ielftype[ie]]=='P2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
        posep = np.where(elftp[ielftype[ie]]=='P3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.0))
        posep = np.where(elftp[ielftype[ie]]=='P4')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.0))
        posep = np.where(elftp[ielftype[ie]]=='P5')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
        ename = 'E2C2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en4PR')
        ielftype = arrset(ielftype,ie,iet_["en4PR"])
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V4')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='P1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
        posep = np.where(elftp[ielftype[ie]]=='P2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.0))
        posep = np.where(elftp[ielftype[ie]]=='P3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.0))
        posep = np.where(elftp[ielftype[ie]]=='P4')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(2.0))
        ename = 'E3C2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en5PR')
        ielftype = arrset(ielftype,ie,iet_["en5PR"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V4')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V5')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='P1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.0))
        posep = np.where(elftp[ielftype[ie]]=='P2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(0.5))
        posep = np.where(elftp[ielftype[ie]]=='P3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-2.0))
        posep = np.where(elftp[ielftype[ie]]=='P4')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.0))
        posep = np.where(elftp[ielftype[ie]]=='P5')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(0.3333333333))
        ename = 'E1C3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en5PR')
        ielftype = arrset(ielftype,ie,iet_["en5PR"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V4')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V5')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='P1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
        posep = np.where(elftp[ielftype[ie]]=='P2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.5))
        posep = np.where(elftp[ielftype[ie]]=='P3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
        posep = np.where(elftp[ielftype[ie]]=='P4')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.0))
        posep = np.where(elftp[ielftype[ie]]=='P5')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(0.3333333333))
        ename = 'E2C3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en5PR')
        ielftype = arrset(ielftype,ie,iet_["en5PR"])
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V4')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V5')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='P1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
        posep = np.where(elftp[ielftype[ie]]=='P2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-0.5))
        posep = np.where(elftp[ielftype[ie]]=='P3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
        posep = np.where(elftp[ielftype[ie]]=='P4')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.0))
        posep = np.where(elftp[ielftype[ie]]=='P5')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-0.5))
        ename = 'E3C3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en4PR')
        ielftype = arrset(ielftype,ie,iet_["en4PR"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V4')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='P1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.0))
        posep = np.where(elftp[ielftype[ie]]=='P2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
        posep = np.where(elftp[ielftype[ie]]=='P3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(0.5))
        posep = np.where(elftp[ielftype[ie]]=='P4')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
        ename = 'E4C3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en5PR')
        ielftype = arrset(ielftype,ie,iet_["en5PR"])
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V4')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V5')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='P1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-2.0))
        posep = np.where(elftp[ielftype[ie]]=='P2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
        posep = np.where(elftp[ielftype[ie]]=='P3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
        posep = np.where(elftp[ielftype[ie]]=='P4')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.0))
        posep = np.where(elftp[ielftype[ie]]=='P5')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
        ename = 'E1C4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en5PR')
        ielftype = arrset(ielftype,ie,iet_["en5PR"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V4')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V5')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='P1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-2.0))
        posep = np.where(elftp[ielftype[ie]]=='P2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
        posep = np.where(elftp[ielftype[ie]]=='P3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.0))
        posep = np.where(elftp[ielftype[ie]]=='P4')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(0.5))
        posep = np.where(elftp[ielftype[ie]]=='P5')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(0.3333333333))
        ename = 'E2C4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en6PR')
        ielftype = arrset(ielftype,ie,iet_["en6PR"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V4')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V5')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V6')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='P1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(0.5))
        posep = np.where(elftp[ielftype[ie]]=='P2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(2.0))
        posep = np.where(elftp[ielftype[ie]]=='P3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
        posep = np.where(elftp[ielftype[ie]]=='P4')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(0.3333333333))
        posep = np.where(elftp[ielftype[ie]]=='P5')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-0.666666666))
        posep = np.where(elftp[ielftype[ie]]=='P6')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(0.25))
        ename = 'E3C4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en5PR')
        ielftype = arrset(ielftype,ie,iet_["en5PR"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V4')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V5')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='P1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-3.0))
        posep = np.where(elftp[ielftype[ie]]=='P2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-2.0))
        posep = np.where(elftp[ielftype[ie]]=='P3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
        posep = np.where(elftp[ielftype[ie]]=='P4')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
        posep = np.where(elftp[ielftype[ie]]=='P5')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(0.75))
        ename = 'E4C4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en3PR')
        ielftype = arrset(ielftype,ie,iet_["en3PR"])
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='P1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-2.0))
        posep = np.where(elftp[ielftype[ie]]=='P2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
        posep = np.where(elftp[ielftype[ie]]=='P3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(0.5))
        ename = 'E1C5'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en5PR')
        ielftype = arrset(ielftype,ie,iet_["en5PR"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V4')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V5')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='P1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
        posep = np.where(elftp[ielftype[ie]]=='P2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.0))
        posep = np.where(elftp[ielftype[ie]]=='P3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(2.0))
        posep = np.where(elftp[ielftype[ie]]=='P4')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-3.0))
        posep = np.where(elftp[ielftype[ie]]=='P5')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(v_['A101']))
        ename = 'E2C5'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en6PR')
        ielftype = arrset(ielftype,ie,iet_["en6PR"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V4')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V5')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V6')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='P1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.0))
        posep = np.where(elftp[ielftype[ie]]=='P2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-2.0))
        posep = np.where(elftp[ielftype[ie]]=='P3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
        posep = np.where(elftp[ielftype[ie]]=='P4')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
        posep = np.where(elftp[ielftype[ie]]=='P5')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.0))
        posep = np.where(elftp[ielftype[ie]]=='P6')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-0.5))
        ename = 'E3C5'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en5PR')
        ielftype = arrset(ielftype,ie,iet_["en5PR"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V4')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V5')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='P1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-2.0))
        posep = np.where(elftp[ielftype[ie]]=='P2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
        posep = np.where(elftp[ielftype[ie]]=='P3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.0))
        posep = np.where(elftp[ielftype[ie]]=='P4')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-2.0))
        posep = np.where(elftp[ielftype[ie]]=='P5')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
        ename = 'E4C5'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en6PR')
        ielftype = arrset(ielftype,ie,iet_["en6PR"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V4')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V5')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(10.0),float(6.0))
        posev = np.where(elftv[ielftype[ie]]=='V6')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='P1')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(2.0))
        posep = np.where(elftp[ielftype[ie]]=='P2')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(2.0))
        posep = np.where(elftp[ielftype[ie]]=='P3')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.0))
        posep = np.where(elftp[ielftype[ie]]=='P4')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(0.5))
        posep = np.where(elftp[ielftype[ie]]=='P5')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-2.0))
        posep = np.where(elftp[ielftype[ie]]=='P6')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        ig = ig_['OBJ']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E1C5'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(10.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E2C5'])
        self.grelw = loaset(self.grelw,ig,posel,float(15.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E3C5'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(20.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E4C5'])
        self.grelw = loaset(self.grelw,ig,posel,float(25.0))
        ig = ig_['CONSTR1']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E1C1'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.5))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E2C1'])
        self.grelw = loaset(self.grelw,ig,posel,float(0.7))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E3C1'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.2))
        ig = ig_['CONSTR2']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E1C2'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.3))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E2C2'])
        self.grelw = loaset(self.grelw,ig,posel,float(0.8))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E3C2'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(3.1))
        ig = ig_['CONSTR3']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E1C3'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(2.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E2C3'])
        self.grelw = loaset(self.grelw,ig,posel,float(0.1))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E3C3'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E4C3'])
        self.grelw = loaset(self.grelw,ig,posel,float(0.65))
        ig = ig_['CONSTR4']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E1C4'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.2))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E2C4'])
        self.grelw = loaset(self.grelw,ig,posel,float(0.3))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E3C4'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.4))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E4C4'])
        self.grelw = loaset(self.grelw,ig,posel,float(0.5))
        ig = ig_['CONSTR5']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E1C5'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(10.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E2C5'])
        self.grelw = loaset(self.grelw,ig,posel,float(15.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E3C5'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(20.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E4C5'])
        self.grelw = loaset(self.grelw,ig,posel,float(25.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               1809.76476
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle)] = grange[legrps]
        self.cupper[np.arange(self.nle)] = np.zeros((self.nle,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-COOR2-AN-7-5"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def en3PR(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        FVALUE  = (
              (EV_[0,0]**self.elpar[iel_][0])*(EV_[1,0]**self.elpar[iel_][1])*(EV_[2,0]**self.elpar[iel_][2]))
        f_   = FVALUE
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = FVALUE*(self.elpar[iel_][0]/EV_[0,0])
            g_[1] = FVALUE*(self.elpar[iel_][1]/EV_[1,0])
            g_[2] = FVALUE*(self.elpar[iel_][2]/EV_[2,0])
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0]  = (
                      FVALUE*(self.elpar[iel_][0]/EV_[0,0])*((self.elpar[iel_][0]-1.0)/EV_[0,0]))
                H_[1,1]  = (
                      FVALUE*(self.elpar[iel_][1]/EV_[1,0])*((self.elpar[iel_][1]-1.0)/EV_[1,0]))
                H_[2,2]  = (
                      FVALUE*(self.elpar[iel_][2]/EV_[2,0])*((self.elpar[iel_][2]-1.0)/EV_[2,0]))
                H_[0,1]  = (
                      FVALUE*(self.elpar[iel_][0]/EV_[0,0])*(self.elpar[iel_][1]/EV_[1,0]))
                H_[1,0] = H_[0,1]
                H_[0,2]  = (
                      FVALUE*(self.elpar[iel_][0]/EV_[0,0])*(self.elpar[iel_][2]/EV_[2,0]))
                H_[2,0] = H_[0,2]
                H_[1,2]  = (
                      FVALUE*(self.elpar[iel_][1]/EV_[1,0])*(self.elpar[iel_][2]/EV_[2,0]))
                H_[2,1] = H_[1,2]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def en4PR(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        FVALUE  = (
              (EV_[0,0]**self.elpar[iel_][0])*(EV_[1,0]**self.elpar[iel_][1])*(EV_[2,0]**self.elpar[iel_][2])*(EV_[3,0]**self.elpar[iel_][3]))
        f_   = FVALUE
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = FVALUE*(self.elpar[iel_][0]/EV_[0,0])
            g_[1] = FVALUE*(self.elpar[iel_][1]/EV_[1,0])
            g_[2] = FVALUE*(self.elpar[iel_][2]/EV_[2,0])
            g_[3] = FVALUE*(self.elpar[iel_][3]/EV_[3,0])
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,0]  = (
                      FVALUE*(self.elpar[iel_][0]/EV_[0,0])*((self.elpar[iel_][0]-1.0)/EV_[0,0]))
                H_[1,1]  = (
                      FVALUE*(self.elpar[iel_][1]/EV_[1,0])*((self.elpar[iel_][1]-1.0)/EV_[1,0]))
                H_[2,2]  = (
                      FVALUE*(self.elpar[iel_][2]/EV_[2,0])*((self.elpar[iel_][2]-1.0)/EV_[2,0]))
                H_[3,3]  = (
                      FVALUE*(self.elpar[iel_][3]/EV_[3,0])*((self.elpar[iel_][3]-1.0)/EV_[3,0]))
                H_[0,1]  = (
                      FVALUE*(self.elpar[iel_][0]/EV_[0,0])*(self.elpar[iel_][1]/EV_[1,0]))
                H_[1,0] = H_[0,1]
                H_[0,2]  = (
                      FVALUE*(self.elpar[iel_][0]/EV_[0,0])*(self.elpar[iel_][2]/EV_[2,0]))
                H_[2,0] = H_[0,2]
                H_[0,3]  = (
                      FVALUE*(self.elpar[iel_][0]/EV_[0,0])*(self.elpar[iel_][3]/EV_[3,0]))
                H_[3,0] = H_[0,3]
                H_[1,2]  = (
                      FVALUE*(self.elpar[iel_][1]/EV_[1,0])*(self.elpar[iel_][2]/EV_[2,0]))
                H_[2,1] = H_[1,2]
                H_[1,3]  = (
                      FVALUE*(self.elpar[iel_][1]/EV_[1,0])*(self.elpar[iel_][3]/EV_[3,0]))
                H_[3,1] = H_[1,3]
                H_[2,3]  = (
                      FVALUE*(self.elpar[iel_][2]/EV_[2,0])*(self.elpar[iel_][3]/EV_[3,0]))
                H_[3,2] = H_[2,3]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def en5PR(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        FVALUE  = (
              (EV_[0,0]**self.elpar[iel_][0])*(EV_[1,0]**self.elpar[iel_][1])*(EV_[2,0]**self.elpar[iel_][2])*(EV_[3,0]**self.elpar[iel_][3])*(EV_[4,0]**self.elpar[iel_][4]))
        f_   = FVALUE
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = FVALUE*(self.elpar[iel_][0]/EV_[0,0])
            g_[1] = FVALUE*(self.elpar[iel_][1]/EV_[1,0])
            g_[2] = FVALUE*(self.elpar[iel_][2]/EV_[2,0])
            g_[3] = FVALUE*(self.elpar[iel_][3]/EV_[3,0])
            g_[4] = FVALUE*(self.elpar[iel_][4]/EV_[4,0])
            if nargout>2:
                H_ = np.zeros((5,5))
                H_[0,0]  = (
                      FVALUE*(self.elpar[iel_][0]/EV_[0,0])*((self.elpar[iel_][0]-1.0)/EV_[0,0]))
                H_[1,1]  = (
                      FVALUE*(self.elpar[iel_][1]/EV_[1,0])*((self.elpar[iel_][1]-1.0)/EV_[1,0]))
                H_[2,2]  = (
                      FVALUE*(self.elpar[iel_][2]/EV_[2,0])*((self.elpar[iel_][2]-1.0)/EV_[2,0]))
                H_[3,3]  = (
                      FVALUE*(self.elpar[iel_][3]/EV_[3,0])*((self.elpar[iel_][3]-1.0)/EV_[3,0]))
                H_[4,4]  = (
                      FVALUE*(self.elpar[iel_][4]/EV_[4,0])*((self.elpar[iel_][4]-1.0)/EV_[4,0]))
                H_[0,1]  = (
                      FVALUE*(self.elpar[iel_][0]/EV_[0,0])*(self.elpar[iel_][1]/EV_[1,0]))
                H_[1,0] = H_[0,1]
                H_[0,2]  = (
                      FVALUE*(self.elpar[iel_][0]/EV_[0,0])*(self.elpar[iel_][2]/EV_[2,0]))
                H_[2,0] = H_[0,2]
                H_[0,3]  = (
                      FVALUE*(self.elpar[iel_][0]/EV_[0,0])*(self.elpar[iel_][3]/EV_[3,0]))
                H_[3,0] = H_[0,3]
                H_[0,4]  = (
                      FVALUE*(self.elpar[iel_][0]/EV_[0,0])*(self.elpar[iel_][4]/EV_[4,0]))
                H_[4,0] = H_[0,4]
                H_[1,2]  = (
                      FVALUE*(self.elpar[iel_][1]/EV_[1,0])*(self.elpar[iel_][2]/EV_[2,0]))
                H_[2,1] = H_[1,2]
                H_[1,3]  = (
                      FVALUE*(self.elpar[iel_][1]/EV_[1,0])*(self.elpar[iel_][3]/EV_[3,0]))
                H_[3,1] = H_[1,3]
                H_[1,4]  = (
                      FVALUE*(self.elpar[iel_][1]/EV_[1,0])*(self.elpar[iel_][4]/EV_[4,0]))
                H_[4,1] = H_[1,4]
                H_[2,3]  = (
                      FVALUE*(self.elpar[iel_][2]/EV_[2,0])*(self.elpar[iel_][3]/EV_[3,0]))
                H_[3,2] = H_[2,3]
                H_[2,4]  = (
                      FVALUE*(self.elpar[iel_][2]/EV_[2,0])*(self.elpar[iel_][4]/EV_[4,0]))
                H_[4,2] = H_[2,4]
                H_[3,4]  = (
                      FVALUE*(self.elpar[iel_][3]/EV_[3,0])*(self.elpar[iel_][4]/EV_[4,0]))
                H_[4,3] = H_[3,4]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def en6PR(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        FVALUE  = (
              (EV_[0,0]**self.elpar[iel_][0])*(EV_[1,0]**self.elpar[iel_][2])*(EV_[2,0]**self.elpar[iel_][3])*(EV_[3,0]**self.elpar[iel_][4])*(EV_[4,0]**self.elpar[iel_][5])*(EV_[5,0]**self.elpar[iel_][1]))
        f_   = FVALUE
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = FVALUE*(self.elpar[iel_][0]/EV_[0,0])
            g_[1] = FVALUE*(self.elpar[iel_][2]/EV_[1,0])
            g_[2] = FVALUE*(self.elpar[iel_][3]/EV_[2,0])
            g_[3] = FVALUE*(self.elpar[iel_][4]/EV_[3,0])
            g_[4] = FVALUE*(self.elpar[iel_][5]/EV_[4,0])
            g_[5] = FVALUE*(self.elpar[iel_][1]/EV_[5,0])
            if nargout>2:
                H_ = np.zeros((6,6))
                H_[0,0]  = (
                      FVALUE*(self.elpar[iel_][0]/EV_[0,0])*((self.elpar[iel_][0]-1.0)/EV_[0,0]))
                H_[1,1]  = (
                      FVALUE*(self.elpar[iel_][2]/EV_[1,0])*((self.elpar[iel_][2]-1.0)/EV_[1,0]))
                H_[2,2]  = (
                      FVALUE*(self.elpar[iel_][3]/EV_[2,0])*((self.elpar[iel_][3]-1.0)/EV_[2,0]))
                H_[3,3]  = (
                      FVALUE*(self.elpar[iel_][4]/EV_[3,0])*((self.elpar[iel_][4]-1.0)/EV_[3,0]))
                H_[4,4]  = (
                      FVALUE*(self.elpar[iel_][5]/EV_[4,0])*((self.elpar[iel_][5]-1.0)/EV_[4,0]))
                H_[5,5]  = (
                      FVALUE*(self.elpar[iel_][1]/EV_[5,0])*((self.elpar[iel_][1]-1.0)/EV_[5,0]))
                H_[0,1]  = (
                      FVALUE*(self.elpar[iel_][0]/EV_[0,0])*(self.elpar[iel_][2]/EV_[1,0]))
                H_[1,0] = H_[0,1]
                H_[0,2]  = (
                      FVALUE*(self.elpar[iel_][0]/EV_[0,0])*(self.elpar[iel_][3]/EV_[2,0]))
                H_[2,0] = H_[0,2]
                H_[0,3]  = (
                      FVALUE*(self.elpar[iel_][0]/EV_[0,0])*(self.elpar[iel_][4]/EV_[3,0]))
                H_[3,0] = H_[0,3]
                H_[0,4]  = (
                      FVALUE*(self.elpar[iel_][0]/EV_[0,0])*(self.elpar[iel_][5]/EV_[4,0]))
                H_[4,0] = H_[0,4]
                H_[0,5]  = (
                      FVALUE*(self.elpar[iel_][0]/EV_[0,0])*(self.elpar[iel_][1]/EV_[5,0]))
                H_[5,0] = H_[0,5]
                H_[1,2]  = (
                      FVALUE*(self.elpar[iel_][2]/EV_[1,0])*(self.elpar[iel_][3]/EV_[2,0]))
                H_[2,1] = H_[1,2]
                H_[1,3]  = (
                      FVALUE*(self.elpar[iel_][2]/EV_[1,0])*(self.elpar[iel_][4]/EV_[3,0]))
                H_[3,1] = H_[1,3]
                H_[1,4]  = (
                      FVALUE*(self.elpar[iel_][2]/EV_[1,0])*(self.elpar[iel_][5]/EV_[4,0]))
                H_[4,1] = H_[1,4]
                H_[1,5]  = (
                      FVALUE*(self.elpar[iel_][2]/EV_[1,0])*(self.elpar[iel_][1]/EV_[5,0]))
                H_[5,1] = H_[1,5]
                H_[2,3]  = (
                      FVALUE*(self.elpar[iel_][3]/EV_[2,0])*(self.elpar[iel_][4]/EV_[3,0]))
                H_[3,2] = H_[2,3]
                H_[2,4]  = (
                      FVALUE*(self.elpar[iel_][3]/EV_[2,0])*(self.elpar[iel_][5]/EV_[4,0]))
                H_[4,2] = H_[2,4]
                H_[2,5]  = (
                      FVALUE*(self.elpar[iel_][3]/EV_[2,0])*(self.elpar[iel_][1]/EV_[5,0]))
                H_[5,2] = H_[2,5]
                H_[3,4]  = (
                      FVALUE*(self.elpar[iel_][4]/EV_[3,0])*(self.elpar[iel_][5]/EV_[4,0]))
                H_[4,3] = H_[3,4]
                H_[3,5]  = (
                      FVALUE*(self.elpar[iel_][4]/EV_[3,0])*(self.elpar[iel_][1]/EV_[5,0]))
                H_[5,3] = H_[3,5]
                H_[4,5]  = (
                      FVALUE*(self.elpar[iel_][5]/EV_[4,0])*(self.elpar[iel_][1]/EV_[5,0]))
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

