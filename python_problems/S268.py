from s2mpjlib import *
class  S268(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : S268
#    *********
# 
#    A quadratic programming problem.
# 
#    Source:
#    K. Schittkowski
#    "More Test Examples for Nonlinear Programming Codes"
#    Springer Verlag, Berlin, Lecture notes in economics and 
#    mathematical systems, volume 282, 1987
# 
#    SIF input: Michel Bierlaire and Annick Sartenaer,
#    October 1992.
# 
#    classification = "QLR2-AN-5-5"
# 
#   the number of functions
# 
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'S268'

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
        v_['D1,1'] = 10197.0
        v_['D1,2'] = -12454.0
        v_['D1,3'] = -1013.0
        v_['D1,4'] = 1948.0
        v_['D1,5'] = 329.0
        v_['D2,1'] = -12454.0
        v_['D2,2'] = 20909.0
        v_['D2,3'] = -1733.0
        v_['D2,4'] = -4914.0
        v_['D2,5'] = -186.0
        v_['D3,1'] = -1013.0
        v_['D3,2'] = -1733.0
        v_['D3,3'] = 1755.0
        v_['D3,4'] = 1089.0
        v_['D3,5'] = -174.0
        v_['D4,1'] = 1948.0
        v_['D4,2'] = -4914.0
        v_['D4,3'] = 1089.0
        v_['D4,4'] = 1515.0
        v_['D4,5'] = -22.0
        v_['D5,1'] = 329.0
        v_['D5,2'] = -186.0
        v_['D5,3'] = -174.0
        v_['D5,4'] = -22.0
        v_['D5,5'] = 27.0
        v_['B1'] = -9170
        v_['B2'] = 17099
        v_['B3'] = -2271
        v_['B4'] = -4336
        v_['B5'] = -43
        v_['1'] = 1
        v_['5'] = 5
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['5'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('NONL',ig_)
        gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['1']),int(v_['5'])+1):
            [ig,ig_,_] = s2mpj_ii('LINEAR',ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = float(v_['B'+str(I)])+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('LINEAR',ig_)
        gtype = arrset(gtype,ig,'<>')
        pbm.gscale = arrset(pbm.gscale,ig,float(-0.5))
        [ig,ig_,_] = s2mpj_ii('C1',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C1')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X2']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X3']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X4']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X5']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C2',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C2')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(10.0)+pbm.A[ig,iv]
        iv = ix_['X2']
        pbm.A[ig,iv] = float(10.0)+pbm.A[ig,iv]
        iv = ix_['X3']
        pbm.A[ig,iv] = float(-3.0)+pbm.A[ig,iv]
        iv = ix_['X4']
        pbm.A[ig,iv] = float(5.0)+pbm.A[ig,iv]
        iv = ix_['X5']
        pbm.A[ig,iv] = float(4.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C3',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C3')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(-8.0)+pbm.A[ig,iv]
        iv = ix_['X2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X3']
        pbm.A[ig,iv] = float(-2.0)+pbm.A[ig,iv]
        iv = ix_['X4']
        pbm.A[ig,iv] = float(-5.0)+pbm.A[ig,iv]
        iv = ix_['X5']
        pbm.A[ig,iv] = float(3.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C4',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C4')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(8.0)+pbm.A[ig,iv]
        iv = ix_['X2']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X3']
        pbm.A[ig,iv] = float(2.0)+pbm.A[ig,iv]
        iv = ix_['X4']
        pbm.A[ig,iv] = float(5.0)+pbm.A[ig,iv]
        iv = ix_['X5']
        pbm.A[ig,iv] = float(-3.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C5',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C5')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(-4.0)+pbm.A[ig,iv]
        iv = ix_['X2']
        pbm.A[ig,iv] = float(-2.0)+pbm.A[ig,iv]
        iv = ix_['X3']
        pbm.A[ig,iv] = float(3.0)+pbm.A[ig,iv]
        iv = ix_['X4']
        pbm.A[ig,iv] = float(-5.0)+pbm.A[ig,iv]
        iv = ix_['X5']
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
        pbm.gconst = arrset(pbm.gconst,ig_['C1'],float(-5.0))
        pbm.gconst = arrset(pbm.gconst,ig_['C2'],float(20.0))
        pbm.gconst = arrset(pbm.gconst,ig_['C3'],float(-40.0))
        pbm.gconst = arrset(pbm.gconst,ig_['C4'],float(11.0))
        pbm.gconst = arrset(pbm.gconst,ig_['C5'],float(-30.0))
        pbm.gconst = arrset(pbm.gconst,ig_['NONL'],float(-14463.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower = np.zeros((pb.n,1))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(1.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'ePROD', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        elftp = []
        elftp = loaset(elftp,it,0,'D')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['1']),int(v_['5'])+1):
            for J in range(int(v_['1']),int(v_['5'])+1):
                ename = 'E'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'ePROD')
                ielftype = arrset(ielftype, ie, iet_["ePROD"])
                vname = 'X'+str(I)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,1.0)
                posev = find(elftv[ielftype[ie]],lambda x:x=='X')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,1.0)
                posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                posep = find(elftp[ielftype[ie]],lambda x:x=='D')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['D'+str(I)+','+str(J)]))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['5'])+1):
            for J in range(int(v_['1']),int(v_['5'])+1):
                ig = ig_['NONL']
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle+pb.neq,pb.m)] = np.zeros((pb.nge,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "QLR2-AN-5-5"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def ePROD(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = pbm.elpar[iel_][0]*EV_[0]*EV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = pbm.elpar[iel_][0]*EV_[1]
            g_[1] = pbm.elpar[iel_][0]*EV_[0]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 0.0
                H_[1,1] = 0.0
                H_[0,1] = pbm.elpar[iel_][0]
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

