from s2xlib import *
class  HS68(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS68
#    *********
# 
#    This is a cost optimal inspection plan.
# 
#    Source: problem 68 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: Nick Gould, August 1991.
# 
#    classification = "OOR2-MN-4-2"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS68'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'HS68'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 4
        v_['A'] = 0.0001
        v_['B'] = 1.0
        v_['D'] = 1.0
        v_['NN'] = 24.0
        v_['1'] = 1
        v_['AN'] = v_['A']*v_['NN']
        v_['ROOTN'] = np.sqrt(v_['NN'])
        v_['DROOTN'] = v_['D']*v_['ROOTN']
        v_['-DROOTN'] = -1.0*v_['DROOTN']
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
        [ig,ig_,_] = s2x_ii('C1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C1')
        iv = ix_['X3']
        pbm.A[ig,iv] = 1.0e+0+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('C2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C2')
        iv = ix_['X4']
        pbm.A[ig,iv] = 1.0e+0+pbm.A[ig,iv]
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
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('inf'))
        pb.xupper = np.full((pb.n,1),float('inf'))
        pb.xlower[ix_['X1']] = 0.0001
        pb.xupper[ix_['X1']] = 100.0
        pb.xupper[ix_['X2']] = 100.0
        pb.xupper[ix_['X3']] = 2.0
        pb.xupper[ix_['X4']] = 2.0
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),1.0)
        pb.y0 = np.full((pb.m,1),1.0)
        pass
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'RECIP', iet_)
        elftv = loaset(elftv,it,0,'X1')
        [it,iet_,_] = s2x_ii( 'NASTYEXP', iet_)
        elftv = loaset(elftv,it,0,'X1')
        elftv = loaset(elftv,it,1,'X3')
        elftv = loaset(elftv,it,2,'X4')
        elftp = []
        elftp = loaset(elftp,it,0,'B')
        [it,iet_,_] = s2x_ii( 'PHI', iet_)
        elftv = loaset(elftv,it,0,'X2')
        elftp = loaset(elftp,it,0,'P')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        ename = 'OE1'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'RECIP')
        ielftype = arrset(ielftype, ie, iet_["RECIP"])
        vname = 'X1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'OE2'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'NASTYEXP')
        ielftype = arrset(ielftype, ie, iet_["NASTYEXP"])
        vname = 'X1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X4')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='B')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],v_['B'])
        ename = 'C1E1'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'PHI')
        ielftype = arrset(ielftype, ie, iet_["PHI"])
        vname = 'X2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],0.0e+0)
        ename = 'C2E1'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'PHI')
        ielftype = arrset(ielftype, ie, iet_["PHI"])
        vname = 'X2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],v_['DROOTN'])
        ename = 'C2E2'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'PHI')
        ielftype = arrset(ielftype, ie, iet_["PHI"])
        vname = 'X2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],v_['-DROOTN'])
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
        pbm.grelw = loaset(pbm.grelw,ig,posel,v_['AN'])
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['OE2'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,-1.0e+0)
        ig = ig_['C1']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C1E1'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,-2.0e+0)
        ig = ig_['C2']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C2E1'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,-1.0e+0)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C2E2'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,-1.0e+0)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "OOR2-MN-4-2"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%
