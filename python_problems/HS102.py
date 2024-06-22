from s2mpjlib import *
class  HS102(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    Source: problem 102 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: N. Gould, December 1989.
# 
#    classification = "OOR2-AN-7-5"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS102'

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
        v_['M'] = 5
        v_['N'] = 7
        v_['A101'] = 0.125
        v_['A102'] = 0.125
        v_['A103'] = 0.5
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        pb.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['1']),int(v_['M'])+1):
            [ig,ig_,_] = s2mpj_ii('CONSTR'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'CONSTR'+str(I))
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
        pbm.gconst = arrset(pbm.gconst,ig_['CONSTR1'],float(1.0))
        pbm.gconst = arrset(pbm.gconst,ig_['CONSTR2'],float(1.0))
        pbm.gconst = arrset(pbm.gconst,ig_['CONSTR3'],float(1.0))
        pbm.gconst = arrset(pbm.gconst,ig_['CONSTR4'],float(1.0))
        pbm.gconst = arrset(pbm.gconst,ig_['CONSTR5'],float(3000.0))
        #%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange = np.full((ngrp,1),None)
        grange[legrps] = np.full((pb.nle,1),float('inf'))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),0.1)
        pb.xupper = np.full((pb.n,1),10.0)
        pb.xlower[ix_['X7']] = 0.01
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(6.0))
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
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        ename = 'E1C1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en4PR')
        ielftype = arrset(ielftype, ie, iet_["en4PR"])
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V4')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.5))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-2.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P4')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        ename = 'E2C1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en5PR')
        ielftype = arrset(ielftype, ie, iet_["en5PR"])
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V4')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V5')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(3.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-2.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P4')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P5')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.5))
        ename = 'E3C1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en5PR')
        ielftype = arrset(ielftype, ie, iet_["en5PR"])
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V4')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V5')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-0.5))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P4')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.66666666))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P5')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.25))
        ename = 'E1C2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en5PR')
        ielftype = arrset(ielftype, ie, iet_["en5PR"])
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V4')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V5')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-0.5))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P4')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P5')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        ename = 'E2C2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en4PR')
        ielftype = arrset(ielftype, ie, iet_["en4PR"])
        vname = 'X3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V4')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P4')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(2.0))
        ename = 'E3C2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en5PR')
        ielftype = arrset(ielftype, ie, iet_["en5PR"])
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V4')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V5')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.5))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-2.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P4')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P5')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.3333333333))
        ename = 'E1C3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en5PR')
        ielftype = arrset(ielftype, ie, iet_["en5PR"])
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V4')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V5')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.5))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P4')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P5')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.3333333333))
        ename = 'E2C3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en5PR')
        ielftype = arrset(ielftype, ie, iet_["en5PR"])
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V4')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V5')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-0.5))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P4')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P5')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-0.5))
        ename = 'E3C3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en4PR')
        ielftype = arrset(ielftype, ie, iet_["en4PR"])
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V4')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.5))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P4')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        ename = 'E4C3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en5PR')
        ielftype = arrset(ielftype, ie, iet_["en5PR"])
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V4')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V5')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-2.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P4')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P5')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        ename = 'E1C4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en5PR')
        ielftype = arrset(ielftype, ie, iet_["en5PR"])
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V4')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V5')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-2.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P4')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.5))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P5')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.3333333333))
        ename = 'E2C4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en6PR')
        ielftype = arrset(ielftype, ie, iet_["en6PR"])
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V4')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V5')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V6')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.5))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(2.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P4')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.3333333333))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P5')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-0.666666666))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P6')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.25))
        ename = 'E3C4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en5PR')
        ielftype = arrset(ielftype, ie, iet_["en5PR"])
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V4')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V5')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-3.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-2.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P4')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P5')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.75))
        ename = 'E4C4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en3PR')
        ielftype = arrset(ielftype, ie, iet_["en3PR"])
        vname = 'X3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-2.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.5))
        ename = 'E1C5'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en5PR')
        ielftype = arrset(ielftype, ie, iet_["en5PR"])
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V4')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V5')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(2.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P4')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-3.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P5')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['A101']))
        ename = 'E2C5'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en6PR')
        ielftype = arrset(ielftype, ie, iet_["en6PR"])
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V4')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V5')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V6')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-2.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P4')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P5')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P6')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-0.5))
        ename = 'E3C5'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en5PR')
        ielftype = arrset(ielftype, ie, iet_["en5PR"])
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V4')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V5')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-2.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P4')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-2.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P5')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        ename = 'E4C5'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en6PR')
        ielftype = arrset(ielftype, ie, iet_["en6PR"])
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V4')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V5')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.1,10.0,6.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V6')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(2.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(2.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P4')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.5))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P5')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-2.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P6')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        ig = ig_['OBJ']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E1C5'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(10.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E2C5'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(15.0))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E3C5'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(20.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E4C5'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(25.0))
        ig = ig_['CONSTR1']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E1C1'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.5))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E2C1'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.7))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E3C1'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.2))
        ig = ig_['CONSTR2']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E1C2'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.3))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E2C2'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.8))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E3C2'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(3.1))
        ig = ig_['CONSTR3']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E1C3'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(2.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E2C3'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.1))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E3C3'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E4C3'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.65))
        ig = ig_['CONSTR4']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E1C4'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.2))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E2C4'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.3))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E3C4'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.4))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E4C4'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.5))
        ig = ig_['CONSTR5']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E1C5'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(10.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E2C5'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(15.0))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E3C5'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(20.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E4C5'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(25.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               1809.76476
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle)] = grange[legrps]
        pb.cupper[np.arange(pb.nle)] = np.zeros((pb.nle,1))
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "OOR2-AN-7-5"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def en3PR(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        FVALUE  = (
              (EV_[0]**pbm.elpar[iel_][0])*(EV_[1]**pbm.elpar[iel_][1])*(EV_[2]**pbm.elpar[iel_][2]))
        f_   = FVALUE
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = FVALUE*(pbm.elpar[iel_][0]/EV_[0])
            g_[1] = FVALUE*(pbm.elpar[iel_][1]/EV_[1])
            g_[2] = FVALUE*(pbm.elpar[iel_][2]/EV_[2])
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0]  = (
                      FVALUE*(pbm.elpar[iel_][0]/EV_[0])*((pbm.elpar[iel_][0]-1.0)/EV_[0]))
                H_[1,1]  = (
                      FVALUE*(pbm.elpar[iel_][1]/EV_[1])*((pbm.elpar[iel_][1]-1.0)/EV_[1]))
                H_[2,2]  = (
                      FVALUE*(pbm.elpar[iel_][2]/EV_[2])*((pbm.elpar[iel_][2]-1.0)/EV_[2]))
                H_[0,1] = FVALUE*(pbm.elpar[iel_][0]/EV_[0])*(pbm.elpar[iel_][1]/EV_[1])
                H_[1,0] = H_[0,1]
                H_[0,2] = FVALUE*(pbm.elpar[iel_][0]/EV_[0])*(pbm.elpar[iel_][2]/EV_[2])
                H_[2,0] = H_[0,2]
                H_[1,2] = FVALUE*(pbm.elpar[iel_][1]/EV_[1])*(pbm.elpar[iel_][2]/EV_[2])
                H_[2,1] = H_[1,2]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def en4PR(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        FVALUE  = (
              (EV_[0]**pbm.elpar[iel_][0])*(EV_[1]**pbm.elpar[iel_][1])*(EV_[2]**pbm.elpar[iel_][2])*(EV_[3]**pbm.elpar[iel_][3]))
        f_   = FVALUE
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = FVALUE*(pbm.elpar[iel_][0]/EV_[0])
            g_[1] = FVALUE*(pbm.elpar[iel_][1]/EV_[1])
            g_[2] = FVALUE*(pbm.elpar[iel_][2]/EV_[2])
            g_[3] = FVALUE*(pbm.elpar[iel_][3]/EV_[3])
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,0]  = (
                      FVALUE*(pbm.elpar[iel_][0]/EV_[0])*((pbm.elpar[iel_][0]-1.0)/EV_[0]))
                H_[1,1]  = (
                      FVALUE*(pbm.elpar[iel_][1]/EV_[1])*((pbm.elpar[iel_][1]-1.0)/EV_[1]))
                H_[2,2]  = (
                      FVALUE*(pbm.elpar[iel_][2]/EV_[2])*((pbm.elpar[iel_][2]-1.0)/EV_[2]))
                H_[3,3]  = (
                      FVALUE*(pbm.elpar[iel_][3]/EV_[3])*((pbm.elpar[iel_][3]-1.0)/EV_[3]))
                H_[0,1] = FVALUE*(pbm.elpar[iel_][0]/EV_[0])*(pbm.elpar[iel_][1]/EV_[1])
                H_[1,0] = H_[0,1]
                H_[0,2] = FVALUE*(pbm.elpar[iel_][0]/EV_[0])*(pbm.elpar[iel_][2]/EV_[2])
                H_[2,0] = H_[0,2]
                H_[0,3] = FVALUE*(pbm.elpar[iel_][0]/EV_[0])*(pbm.elpar[iel_][3]/EV_[3])
                H_[3,0] = H_[0,3]
                H_[1,2] = FVALUE*(pbm.elpar[iel_][1]/EV_[1])*(pbm.elpar[iel_][2]/EV_[2])
                H_[2,1] = H_[1,2]
                H_[1,3] = FVALUE*(pbm.elpar[iel_][1]/EV_[1])*(pbm.elpar[iel_][3]/EV_[3])
                H_[3,1] = H_[1,3]
                H_[2,3] = FVALUE*(pbm.elpar[iel_][2]/EV_[2])*(pbm.elpar[iel_][3]/EV_[3])
                H_[3,2] = H_[2,3]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def en5PR(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        FVALUE  = (
              (EV_[0]**pbm.elpar[iel_][0])*(EV_[1]**pbm.elpar[iel_][1])*(EV_[2]**pbm.elpar[iel_][2])*(EV_[3]**pbm.elpar[iel_][3])*(EV_[4]**pbm.elpar[iel_][4]))
        f_   = FVALUE
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = FVALUE*(pbm.elpar[iel_][0]/EV_[0])
            g_[1] = FVALUE*(pbm.elpar[iel_][1]/EV_[1])
            g_[2] = FVALUE*(pbm.elpar[iel_][2]/EV_[2])
            g_[3] = FVALUE*(pbm.elpar[iel_][3]/EV_[3])
            g_[4] = FVALUE*(pbm.elpar[iel_][4]/EV_[4])
            if nargout>2:
                H_ = np.zeros((5,5))
                H_[0,0]  = (
                      FVALUE*(pbm.elpar[iel_][0]/EV_[0])*((pbm.elpar[iel_][0]-1.0)/EV_[0]))
                H_[1,1]  = (
                      FVALUE*(pbm.elpar[iel_][1]/EV_[1])*((pbm.elpar[iel_][1]-1.0)/EV_[1]))
                H_[2,2]  = (
                      FVALUE*(pbm.elpar[iel_][2]/EV_[2])*((pbm.elpar[iel_][2]-1.0)/EV_[2]))
                H_[3,3]  = (
                      FVALUE*(pbm.elpar[iel_][3]/EV_[3])*((pbm.elpar[iel_][3]-1.0)/EV_[3]))
                H_[4,4]  = (
                      FVALUE*(pbm.elpar[iel_][4]/EV_[4])*((pbm.elpar[iel_][4]-1.0)/EV_[4]))
                H_[0,1] = FVALUE*(pbm.elpar[iel_][0]/EV_[0])*(pbm.elpar[iel_][1]/EV_[1])
                H_[1,0] = H_[0,1]
                H_[0,2] = FVALUE*(pbm.elpar[iel_][0]/EV_[0])*(pbm.elpar[iel_][2]/EV_[2])
                H_[2,0] = H_[0,2]
                H_[0,3] = FVALUE*(pbm.elpar[iel_][0]/EV_[0])*(pbm.elpar[iel_][3]/EV_[3])
                H_[3,0] = H_[0,3]
                H_[0,4] = FVALUE*(pbm.elpar[iel_][0]/EV_[0])*(pbm.elpar[iel_][4]/EV_[4])
                H_[4,0] = H_[0,4]
                H_[1,2] = FVALUE*(pbm.elpar[iel_][1]/EV_[1])*(pbm.elpar[iel_][2]/EV_[2])
                H_[2,1] = H_[1,2]
                H_[1,3] = FVALUE*(pbm.elpar[iel_][1]/EV_[1])*(pbm.elpar[iel_][3]/EV_[3])
                H_[3,1] = H_[1,3]
                H_[1,4] = FVALUE*(pbm.elpar[iel_][1]/EV_[1])*(pbm.elpar[iel_][4]/EV_[4])
                H_[4,1] = H_[1,4]
                H_[2,3] = FVALUE*(pbm.elpar[iel_][2]/EV_[2])*(pbm.elpar[iel_][3]/EV_[3])
                H_[3,2] = H_[2,3]
                H_[2,4] = FVALUE*(pbm.elpar[iel_][2]/EV_[2])*(pbm.elpar[iel_][4]/EV_[4])
                H_[4,2] = H_[2,4]
                H_[3,4] = FVALUE*(pbm.elpar[iel_][3]/EV_[3])*(pbm.elpar[iel_][4]/EV_[4])
                H_[4,3] = H_[3,4]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def en6PR(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        FVALUE  = (
              (EV_[0]**pbm.elpar[iel_][0])*(EV_[1]**pbm.elpar[iel_][2])*(EV_[2]**pbm.elpar[iel_][3])*(EV_[3]**pbm.elpar[iel_][4])*(EV_[4]**pbm.elpar[iel_][5])*(EV_[5]**pbm.elpar[iel_][1]))
        f_   = FVALUE
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = FVALUE*(pbm.elpar[iel_][0]/EV_[0])
            g_[1] = FVALUE*(pbm.elpar[iel_][2]/EV_[1])
            g_[2] = FVALUE*(pbm.elpar[iel_][3]/EV_[2])
            g_[3] = FVALUE*(pbm.elpar[iel_][4]/EV_[3])
            g_[4] = FVALUE*(pbm.elpar[iel_][5]/EV_[4])
            g_[5] = FVALUE*(pbm.elpar[iel_][1]/EV_[5])
            if nargout>2:
                H_ = np.zeros((6,6))
                H_[0,0]  = (
                      FVALUE*(pbm.elpar[iel_][0]/EV_[0])*((pbm.elpar[iel_][0]-1.0)/EV_[0]))
                H_[1,1]  = (
                      FVALUE*(pbm.elpar[iel_][2]/EV_[1])*((pbm.elpar[iel_][2]-1.0)/EV_[1]))
                H_[2,2]  = (
                      FVALUE*(pbm.elpar[iel_][3]/EV_[2])*((pbm.elpar[iel_][3]-1.0)/EV_[2]))
                H_[3,3]  = (
                      FVALUE*(pbm.elpar[iel_][4]/EV_[3])*((pbm.elpar[iel_][4]-1.0)/EV_[3]))
                H_[4,4]  = (
                      FVALUE*(pbm.elpar[iel_][5]/EV_[4])*((pbm.elpar[iel_][5]-1.0)/EV_[4]))
                H_[5,5]  = (
                      FVALUE*(pbm.elpar[iel_][1]/EV_[5])*((pbm.elpar[iel_][1]-1.0)/EV_[5]))
                H_[0,1] = FVALUE*(pbm.elpar[iel_][0]/EV_[0])*(pbm.elpar[iel_][2]/EV_[1])
                H_[1,0] = H_[0,1]
                H_[0,2] = FVALUE*(pbm.elpar[iel_][0]/EV_[0])*(pbm.elpar[iel_][3]/EV_[2])
                H_[2,0] = H_[0,2]
                H_[0,3] = FVALUE*(pbm.elpar[iel_][0]/EV_[0])*(pbm.elpar[iel_][4]/EV_[3])
                H_[3,0] = H_[0,3]
                H_[0,4] = FVALUE*(pbm.elpar[iel_][0]/EV_[0])*(pbm.elpar[iel_][5]/EV_[4])
                H_[4,0] = H_[0,4]
                H_[0,5] = FVALUE*(pbm.elpar[iel_][0]/EV_[0])*(pbm.elpar[iel_][1]/EV_[5])
                H_[5,0] = H_[0,5]
                H_[1,2] = FVALUE*(pbm.elpar[iel_][2]/EV_[1])*(pbm.elpar[iel_][3]/EV_[2])
                H_[2,1] = H_[1,2]
                H_[1,3] = FVALUE*(pbm.elpar[iel_][2]/EV_[1])*(pbm.elpar[iel_][4]/EV_[3])
                H_[3,1] = H_[1,3]
                H_[1,4] = FVALUE*(pbm.elpar[iel_][2]/EV_[1])*(pbm.elpar[iel_][5]/EV_[4])
                H_[4,1] = H_[1,4]
                H_[1,5] = FVALUE*(pbm.elpar[iel_][2]/EV_[1])*(pbm.elpar[iel_][1]/EV_[5])
                H_[5,1] = H_[1,5]
                H_[2,3] = FVALUE*(pbm.elpar[iel_][3]/EV_[2])*(pbm.elpar[iel_][4]/EV_[3])
                H_[3,2] = H_[2,3]
                H_[2,4] = FVALUE*(pbm.elpar[iel_][3]/EV_[2])*(pbm.elpar[iel_][5]/EV_[4])
                H_[4,2] = H_[2,4]
                H_[2,5] = FVALUE*(pbm.elpar[iel_][3]/EV_[2])*(pbm.elpar[iel_][1]/EV_[5])
                H_[5,2] = H_[2,5]
                H_[3,4] = FVALUE*(pbm.elpar[iel_][4]/EV_[3])*(pbm.elpar[iel_][5]/EV_[4])
                H_[4,3] = H_[3,4]
                H_[3,5] = FVALUE*(pbm.elpar[iel_][4]/EV_[3])*(pbm.elpar[iel_][1]/EV_[5])
                H_[5,3] = H_[3,5]
                H_[4,5] = FVALUE*(pbm.elpar[iel_][5]/EV_[4])*(pbm.elpar[iel_][1]/EV_[5])
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

