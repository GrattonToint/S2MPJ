from s2xlib import *
class  HS104(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS104
#    *********
# 
#    Source: problem 104 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: Nick Gould, August 1991.
# 
#    classification = "OOR2-AN-8-5"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS104'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'HS104'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 8
        v_['1'] = 1
        v_['4'] = 4
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2x_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2x_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(-1.0e+0)+pbm.A[ig,iv]
        iv = ix_['X2']
        pbm.A[ig,iv] = float(-1.0e+0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('C1',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'C1')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(1.0e-1)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('C2',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'C2')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(1.0e-1)+pbm.A[ig,iv]
        iv = ix_['X2']
        pbm.A[ig,iv] = float(1.0e-1)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('C3',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'C3')
        [ig,ig_,_] = s2x_ii('C4',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'C4')
        [ig,ig_,_] = s2x_ii('C5',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C5')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(-1.0e+0)+pbm.A[ig,iv]
        iv = ix_['X2']
        pbm.A[ig,iv] = float(-1.0e+0)+pbm.A[ig,iv]
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
        #%%%%%%%%%%%%%%%%%%  CONSTANTS %%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.full((ngrp,1),1.0e+0)
        pbm.gconst = arrset(pbm.gconst,ig_['OBJ'],float(-1.0e+1))
        pbm.gconst = arrset(pbm.gconst,ig_['C5'],float(-9.0e+0))
        #%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange = np.full((ngrp,1),None)
        grange[legrps] = np.full((pb.nle,1),float('inf'))
        grange[gegrps] = np.full((pb.nge,1),float('inf'))
        grange = arrset(grange,ig_['C5'],float(3.2e+0))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),1.0e-1)
        pb.xupper = np.full((pb.n,1),1.0e+1)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        if('X1' in ix_):
            pb.x0[ix_['X1']] = float(6.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X1']),float(6.0)))
        if('X2' in ix_):
            pb.x0[ix_['X2']] = float(3.0)
        else:
            pb.y0 = arrset(pb.y0,find(pbm.congrps,lambda x:x==ig_['X2']),float(3.0))
        if('X3' in ix_):
            pb.x0[ix_['X3']] = float(0.4)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X3']),float(0.4)))
        if('X4' in ix_):
            pb.x0[ix_['X4']] = float(0.2)
        else:
            pb.y0 = arrset(pb.y0,find(pbm.congrps,lambda x:x==ig_['X4']),float(0.2))
        if('X5' in ix_):
            pb.x0[ix_['X5']] = float(6.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X5']),float(6.0)))
        if('X6' in ix_):
            pb.x0[ix_['X6']] = float(6.0)
        else:
            pb.y0 = arrset(pb.y0,find(pbm.congrps,lambda x:x==ig_['X6']),float(6.0))
        if('X7' in ix_):
            pb.x0[ix_['X7']] = float(1.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X7']),float(1.0)))
        if('X8' in ix_):
            pb.x0[ix_['X8']] = float(0.5)
        else:
            pb.y0 = arrset(pb.y0,find(pbm.congrps,lambda x:x==ig_['X8']),float(0.5))
        pass
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'ePROD', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftp = []
        elftp = loaset(elftp,it,0,'P1')
        elftp = loaset(elftp,it,1,'P2')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        ename = 'OE1'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePROD')
        ielftype = arrset(ielftype, ie, iet_["ePROD"])
        vname = 'X1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,1.0e-1,1.0e+1,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,1.0e-1,1.0e+1,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.67))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-0.67))
        ename = 'OE2'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePROD')
        ielftype = arrset(ielftype, ie, iet_["ePROD"])
        vname = 'X2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,1.0e-1,1.0e+1,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X8'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,1.0e-1,1.0e+1,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.67))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-0.67))
        ename = 'C1E1'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePROD')
        ielftype = arrset(ielftype, ie, iet_["ePROD"])
        vname = 'X5'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,1.0e-1,1.0e+1,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,1.0e-1,1.0e+1,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        ename = 'C2E1'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePROD')
        ielftype = arrset(ielftype, ie, iet_["ePROD"])
        vname = 'X6'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,1.0e-1,1.0e+1,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X8'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,1.0e-1,1.0e+1,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        ename = 'C3E1'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePROD')
        ielftype = arrset(ielftype, ie, iet_["ePROD"])
        vname = 'X3'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,1.0e-1,1.0e+1,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,1.0e-1,1.0e+1,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
        ename = 'C3E2'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePROD')
        ielftype = arrset(ielftype, ie, iet_["ePROD"])
        vname = 'X3'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,1.0e-1,1.0e+1,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,1.0e-1,1.0e+1,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-0.71))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
        ename = 'C3E3'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePROD')
        ielftype = arrset(ielftype, ie, iet_["ePROD"])
        vname = 'X3'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,1.0e-1,1.0e+1,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,1.0e-1,1.0e+1,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.3))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        ename = 'C4E1'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePROD')
        ielftype = arrset(ielftype, ie, iet_["ePROD"])
        vname = 'X4'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,1.0e-1,1.0e+1,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,1.0e-1,1.0e+1,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
        ename = 'C4E2'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePROD')
        ielftype = arrset(ielftype, ie, iet_["ePROD"])
        vname = 'X4'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,1.0e-1,1.0e+1,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,1.0e-1,1.0e+1,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-0.71))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
        ename = 'C4E3'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePROD')
        ielftype = arrset(ielftype, ie, iet_["ePROD"])
        vname = 'X4'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,1.0e-1,1.0e+1,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X8'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,1.0e-1,1.0e+1,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.3))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
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
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['OE1'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(4.0e-1))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['OE2'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(4.0e-1))
        ig = ig_['C1']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C1E1'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(5.88e-2))
        ig = ig_['C2']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C2E1'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(5.88e-2))
        ig = ig_['C3']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C3E1'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(4.0e+0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C3E2'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(2.0e+0))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C3E3'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(5.88e-2))
        ig = ig_['C4']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C4E1'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(4.0e+0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C4E2'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(2.0e+0))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C4E3'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(5.88e-2))
        ig = ig_['C5']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['OE1'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(4.0e-1))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['OE2'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(4.0e-1))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle)] = grange[legrps]
        pb.cupper[np.arange(pb.nle)] = np.zeros((pb.nle,1))
        pb.clower[np.arange(pb.nle+pb.neq,pb.m)] = np.zeros((pb.nge,1))
        pb.cupper[np.arange(pb.nge)] = grange[gegrps]
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "OOR2-AN-8-5"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def ePROD(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = (EV_[0]**pbm.elpar[iel_][0])*(EV_[1]**pbm.elpar[iel_][1])
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0]  = (
                  pbm.elpar[iel_][0]*(EV_[0]**(pbm.elpar[iel_][0]-1.0))*(EV_[1]**pbm.elpar[iel_][1]))
            g_[1]  = (
                  pbm.elpar[iel_][1]*(EV_[0]**pbm.elpar[iel_][0])*(EV_[1]**(pbm.elpar[iel_][1]-1.0)))
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0]  = (
                      pbm.elpar[iel_][0]*(EV_[0]**(pbm.elpar[iel_][0]-2.0))*(pbm.elpar[iel_][0]-1.0)*(EV_[1]**pbm.elpar[iel_][1]))
                H_[0,1]  = (
                      pbm.elpar[iel_][0]*(EV_[0]**(pbm.elpar[iel_][0]-1.0))*pbm.elpar[iel_][1]*(EV_[1]**(pbm.elpar[iel_][1]-1.0)))
                H_[1,0] = H_[0,1]
                H_[1,1]  = (
                      pbm.elpar[iel_][1]*(pbm.elpar[iel_][1]-1.0)*(EV_[0]**pbm.elpar[iel_][0])*(EV_[1]**(pbm.elpar[iel_][1]-2.0)))
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

