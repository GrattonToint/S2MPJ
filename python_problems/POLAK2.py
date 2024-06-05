from s2mpjlib import *
class  POLAK2(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : POLAK2
#    *********
# 
#    A nonlinear minmax problem in ten variables.
# 
#    Source: 
#    E. Polak, D.H. Mayne and J.E. Higgins,
#    "Superlinearly convergent algorithm for min-max problems"
#    JOTA 69, pp. 407-439, 1991.
# 
#    SIF input: Ph. Toint, Nov 1993.
# 
#    classification = "LOR2-AN-11-2"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'POLAK2'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'POLAK2'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['1'] = 1
        v_['10'] = 10
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['10'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
        [iv,ix_,_] = s2mpj_ii('U',ix_)
        pb.xnames=arrset(pb.xnames,iv,'U')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['U']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('F1',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'F1')
        iv = ix_['U']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('F2',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'F2')
        iv = ix_['U']
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
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(0.1))
        if('X1' in ix_):
            pb.x0[ix_['X1']] = float(100.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X1']),float(100.0)))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eEL', iet_)
        elftv = loaset(elftv,it,0,'XX1')
        elftv = loaset(elftv,it,1,'XX2')
        elftv = loaset(elftv,it,2,'XX3')
        elftv = loaset(elftv,it,3,'XX4')
        elftv = loaset(elftv,it,4,'XX5')
        elftv = loaset(elftv,it,5,'XX6')
        elftv = loaset(elftv,it,6,'XX7')
        elftv = loaset(elftv,it,7,'XX8')
        elftv = loaset(elftv,it,8,'XX9')
        elftv = loaset(elftv,it,9,'XX10')
        elftp = []
        elftp = loaset(elftp,it,0,'P')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        ename = 'E1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eEL')
        ielftype = arrset(ielftype, ie, iet_["eEL"])
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.1)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.1)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.1)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.1)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX4')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.1)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX5')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.1)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX6')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.1)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX7')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X8'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.1)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX8')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.1)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX9')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.1)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX10')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(2.0))
        ename = 'E2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eEL')
        ielftype = arrset(ielftype, ie, iet_["eEL"])
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.1)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.1)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.1)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.1)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX4')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.1)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX5')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.1)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX6')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.1)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX7')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X8'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.1)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX8')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.1)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX9')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.1)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX10')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-2.0))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        ig = ig_['F1']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E1'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['F2']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E2'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               54.598146
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.cupper[np.arange(pb.nle)] = np.zeros((pb.nle,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "LOR2-AN-11-2"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eEL(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        A = 1.0e-8*EV_[0]*EV_[0]+(EV_[1]+pbm.elpar[iel_][0])**2
        A = A+EV_[2]*EV_[2]+4.0*EV_[3]*EV_[3]
        A = A+EV_[4]*EV_[4]+EV_[5]*EV_[5]+EV_[6]*EV_[6]
        A = A+EV_[7]*EV_[7]+EV_[8]*EV_[8]+EV_[9]*EV_[9]
        EA = np.exp(A)
        f_   = EA
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0e-8*EV_[0]*EA
            g_[1] = 2.0*(EV_[1]+pbm.elpar[iel_][0])*EA
            g_[2] = 2.0*EV_[2]*EA
            g_[3] = 8.0*EV_[3]*EA
            g_[4] = 2.0*EV_[4]*EA
            g_[5] = 2.0*EV_[5]*EA
            g_[6] = 2.0*EV_[6]*EA
            g_[7] = 2.0*EV_[7]*EA
            g_[8] = 2.0*EV_[8]*EA
            g_[9] = 2.0*EV_[9]*EA
            if nargout>2:
                H_ = np.zeros((10,10))
                H_[0,0] = 2.0e-8*EA*(1.0+2.0e-8*EV_[0]**2)
                H_[0,1] = 4.0e-8*EV_[0]*(EV_[1]+pbm.elpar[iel_][0])*EA
                H_[1,0] = H_[0,1]
                H_[0,2] = 4.0e-8*EV_[0]*EV_[2]*EA
                H_[2,0] = H_[0,2]
                H_[0,3] = 1.6e-7*EV_[0]*EV_[3]*EA
                H_[3,0] = H_[0,3]
                H_[0,4] = 4.0e-8*EV_[0]*EV_[4]*EA
                H_[4,0] = H_[0,4]
                H_[0,5] = 4.0e-8*EV_[0]*EV_[5]*EA
                H_[5,0] = H_[0,5]
                H_[0,6] = 4.0e-8*EV_[0]*EV_[6]*EA
                H_[6,0] = H_[0,6]
                H_[0,7] = 4.0e-8*EV_[0]*EV_[7]*EA
                H_[7,0] = H_[0,7]
                H_[0,8] = 4.0e-8*EV_[0]*EV_[8]*EA
                H_[8,0] = H_[0,8]
                H_[0,9] = 4.0e-8*EV_[0]*EV_[9]*EA
                H_[9,0] = H_[0,9]
                H_[1,1] = 2.0*EA*(1.0+2.0*(EV_[1]+pbm.elpar[iel_][0])**2)
                H_[1,2] = 4.0*(EV_[1]+pbm.elpar[iel_][0])*EV_[2]*EA
                H_[2,1] = H_[1,2]
                H_[1,3] = 16.0*(EV_[1]+pbm.elpar[iel_][0])*EV_[3]*EA
                H_[3,1] = H_[1,3]
                H_[1,4] = 4.0*(EV_[1]+pbm.elpar[iel_][0])*EV_[4]*EA
                H_[4,1] = H_[1,4]
                H_[1,5] = 4.0*(EV_[1]+pbm.elpar[iel_][0])*EV_[5]*EA
                H_[5,1] = H_[1,5]
                H_[1,6] = 4.0*(EV_[1]+pbm.elpar[iel_][0])*EV_[6]*EA
                H_[6,1] = H_[1,6]
                H_[1,7] = 4.0*(EV_[1]+pbm.elpar[iel_][0])*EV_[7]*EA
                H_[7,1] = H_[1,7]
                H_[1,8] = 4.0*(EV_[1]+pbm.elpar[iel_][0])*EV_[8]*EA
                H_[8,1] = H_[1,8]
                H_[1,9] = 4.0*(EV_[1]+pbm.elpar[iel_][0])*EV_[9]*EA
                H_[9,1] = H_[1,9]
                H_[2,2] = 2.0*EA*(1.0+2.0*EV_[2]*EV_[2])
                H_[2,3] = 16.0*EV_[2]*EV_[3]*EA
                H_[3,2] = H_[2,3]
                H_[2,4] = 4.0*EV_[2]*EV_[4]*EA
                H_[4,2] = H_[2,4]
                H_[2,5] = 4.0*EV_[2]*EV_[5]*EA
                H_[5,2] = H_[2,5]
                H_[2,6] = 4.0*EV_[2]*EV_[6]*EA
                H_[6,2] = H_[2,6]
                H_[2,7] = 4.0*EV_[2]*EV_[7]*EA
                H_[7,2] = H_[2,7]
                H_[2,8] = 4.0*EV_[2]*EV_[8]*EA
                H_[8,2] = H_[2,8]
                H_[2,9] = 4.0*EV_[2]*EV_[9]*EA
                H_[9,2] = H_[2,9]
                H_[3,3] = 8.0*EA*(1.0+8.0*EV_[3]*EV_[3])
                H_[3,4] = 16.0*EV_[3]*EV_[4]*EA
                H_[4,3] = H_[3,4]
                H_[3,5] = 16.0*EV_[3]*EV_[5]*EA
                H_[5,3] = H_[3,5]
                H_[3,6] = 16.0*EV_[3]*EV_[6]*EA
                H_[6,3] = H_[3,6]
                H_[3,7] = 16.0*EV_[3]*EV_[7]*EA
                H_[7,3] = H_[3,7]
                H_[3,8] = 16.0*EV_[3]*EV_[8]*EA
                H_[8,3] = H_[3,8]
                H_[3,9] = 16.0*EV_[3]*EV_[9]*EA
                H_[9,3] = H_[3,9]
                H_[4,4] = 2.0*EA*(1.0+2.0*EV_[4]*EV_[4])
                H_[4,5] = 4.0*EV_[4]*EV_[5]*EA
                H_[5,4] = H_[4,5]
                H_[4,6] = 4.0*EV_[4]*EV_[6]*EA
                H_[6,4] = H_[4,6]
                H_[4,7] = 4.0*EV_[4]*EV_[7]*EA
                H_[7,4] = H_[4,7]
                H_[4,8] = 4.0*EV_[4]*EV_[8]*EA
                H_[8,4] = H_[4,8]
                H_[4,9] = 4.0*EV_[4]*EV_[9]*EA
                H_[9,4] = H_[4,9]
                H_[5,5] = 2.0*EA*(1.0+2.0*EV_[5]*EV_[5])
                H_[5,6] = 4.0*EV_[5]*EV_[6]*EA
                H_[6,5] = H_[5,6]
                H_[5,7] = 4.0*EV_[5]*EV_[7]*EA
                H_[7,5] = H_[5,7]
                H_[5,8] = 4.0*EV_[5]*EV_[8]*EA
                H_[8,5] = H_[5,8]
                H_[5,9] = 4.0*EV_[5]*EV_[9]*EA
                H_[9,5] = H_[5,9]
                H_[6,6] = 2.0*EA*(1.0+2.0*EV_[6]*EV_[6])
                H_[6,7] = 4.0*EV_[6]*EV_[7]*EA
                H_[7,6] = H_[6,7]
                H_[6,8] = 4.0*EV_[6]*EV_[8]*EA
                H_[8,6] = H_[6,8]
                H_[6,9] = 4.0*EV_[6]*EV_[9]*EA
                H_[9,6] = H_[6,9]
                H_[7,7] = 2.0*EA*(1.0+2.0*EV_[7]*EV_[7])
                H_[7,8] = 4.0*EV_[7]*EV_[8]*EA
                H_[8,7] = H_[7,8]
                H_[7,9] = 4.0*EV_[7]*EV_[9]*EA
                H_[9,7] = H_[7,9]
                H_[8,8] = 2.0*EA*(1.0+2.0*EV_[8]*EV_[8])
                H_[8,9] = 4.0*EV_[8]*EV_[9]*EA
                H_[9,8] = H_[8,9]
                H_[9,9] = 2.0*EA*(1.0+2.0*EV_[9]*EV_[9])
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

