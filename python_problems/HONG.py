from s2mpjlib import *
class  HONG(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    Source: Se June Hong/Chid Apte
# 
#    SIF input: A.R.Conn, Jan 1991.
# 
#    classification = "OLR2-AN-4-1"
# 
#   Problem parameters
# 
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HONG'

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
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2mpj_ii('T1',ix_)
        pb.xnames=arrset(pb.xnames,iv,'T1')
        [iv,ix_,_] = s2mpj_ii('T2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'T2')
        [iv,ix_,_] = s2mpj_ii('T3',ix_)
        pb.xnames=arrset(pb.xnames,iv,'T3')
        [iv,ix_,_] = s2mpj_ii('T4',ix_)
        pb.xnames=arrset(pb.xnames,iv,'T4')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('SUM1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'SUM1')
        iv = ix_['T1']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['T2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['T3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['T4']
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
        pbm.gconst = arrset(pbm.gconst,ig_['SUM1'],float(1.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),0.0)
        pb.xupper = np.full((pb.n,1),1.0)
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(0.5))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eEXP', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftp = []
        elftp = loaset(elftp,it,0,'P1')
        elftp = loaset(elftp,it,1,'P2')
        elftp = loaset(elftp,it,2,'P3')
        elftp = loaset(elftp,it,3,'P4')
        elftp = loaset(elftp,it,4,'P5')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        ename = 'E1'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eEXP')
            ielftype = arrset( ielftype,ie,iet_['eEXP'])
        vname = 'T1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,0.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
        loaset(pbm.elpar,ie,posep[0],float(25.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.92))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P4')
        loaset(pbm.elpar,ie,posep[0],float(0.08))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P5')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.38))
        ename = 'E2'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eEXP')
            ielftype = arrset( ielftype,ie,iet_['eEXP'])
        vname = 'T2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,0.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
        loaset(pbm.elpar,ie,posep[0],float(50.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-2.95))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P4')
        loaset(pbm.elpar,ie,posep[0],float(3.95))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P5')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.11))
        ename = 'E3'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eEXP')
            ielftype = arrset( ielftype,ie,iet_['eEXP'])
        vname = 'T3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,0.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(9.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
        loaset(pbm.elpar,ie,posep[0],float(-4.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.66))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P4')
        loaset(pbm.elpar,ie,posep[0],float(1657834.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P5')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.48))
        ename = 'E4'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eEXP')
            ielftype = arrset( ielftype,ie,iet_['eEXP'])
        vname = 'T4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,0.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
        loaset(pbm.elpar,ie,posep[0],float(20000.0))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.11))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P4')
        loaset(pbm.elpar,ie,posep[0],float(0.89))
        posep = find(elftp[ielftype[ie]],lambda x:x=='P5')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.00035))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        ig = ig_['OBJ']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E1'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E2'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E3'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E4'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution unknown
        pb.objlower = -4.0
        pb.objupper = 300.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
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
        pb.pbclass = "OLR2-AN-4-1"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eEXP(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        XTOT = pbm.elpar[iel_][0]+pbm.elpar[iel_][1]*EV_[0]
        EP5 = np.exp(pbm.elpar[iel_][4]*XTOT)
        f_   = pbm.elpar[iel_][2]+pbm.elpar[iel_][3]*EP5
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = pbm.elpar[iel_][1]*pbm.elpar[iel_][3]*pbm.elpar[iel_][4]*EP5
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0]  = (
                      pbm.elpar[iel_][1]*pbm.elpar[iel_][1]*pbm.elpar[iel_][3]*pbm.elpar[iel_][4]*pbm.elpar[iel_][4]*EP5)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

