from s2mpjlib import *
class  HIMMELP6(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HIMMELP6
#    *********
# 
#    A nonlinear problem with inequality constraints, attributed to Himmelblau
#    by B.N. Pshenichnyj (case V).
# 
#    The problem is nonconvex and has redundant constraints at the solution.
# 
#    Source: 
#    B.N. Pshenichnyj
#    "The Linearization Method for Constrained Optimization",
#    Springer Verlag, SCM Series 22, Heidelberg, 1994
# 
#    SIF input: Ph. Toint, December 1994.
# 
#    classification = "OQR2-AN-2-5"
# 
#    Problem data
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HIMMELP6'

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
        v_['B1'] = 0.1963666677
        v_['B1'] = 75.0+v_['B1']
        v_['B2'] = -.8112755343
        v_['B2'] = -3.0+v_['B2']
        v_['B6'] = -.8306567613
        v_['B6'] = -6.0+v_['B6']
        v_['-B2'] = -1.0*v_['B2']
        v_['-B6'] = -1.0*v_['B6']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2mpj_ii('X1',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X1')
        [iv,ix_,_] = s2mpj_ii('X2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X2')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(v_['-B2'])+pbm.A[ig,iv]
        iv = ix_['X2']
        pbm.A[ig,iv] = float(v_['-B6'])+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C5',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C5')
        [ig,ig_,_] = s2mpj_ii('C8',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'C8')
        iv = ix_['X2']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C9',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'C9')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(5.0)+pbm.A[ig,iv]
        iv = ix_['X2']
        pbm.A[ig,iv] = float(100.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C11',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'C11')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X2']
        pbm.A[ig,iv] = float(-1.5)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C12',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'C12')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X2']
        pbm.A[ig,iv] = float(0.16)+pbm.A[ig,iv]
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
        pbm.gconst = arrset(pbm.gconst,ig_['OBJ'],float(v_['B1']))
        pbm.gconst = arrset(pbm.gconst,ig_['C5'],float(700.0))
        pbm.gconst = arrset(pbm.gconst,ig_['C9'],float(2775.0))
        pbm.gconst = arrset(pbm.gconst,ig_['C11'],float(22.5))
        pbm.gconst = arrset(pbm.gconst,ig_['C12'],float(41.4))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),float('inf'))
        pb.xupper[ix_['X1']] = 75.0
        pb.xupper[ix_['X2']] = 65.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        pb.x0[ix_['X1']] = float(68.8)
        pb.x0[ix_['X2']] = float(31.2)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eOBNL', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        [it,iet_,_] = s2mpj_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        [it,iet_,_] = s2mpj_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        ename = 'OB'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eOBNL')
        ielftype = arrset(ielftype, ie, iet_["eOBNL"])
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'X1X2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
        ielftype = arrset(ielftype, ie, iet_["en2PR"])
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'X1SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
        ielftype = arrset(ielftype, ie, iet_["eSQ"])
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'X2SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
        ielftype = arrset(ielftype, ie, iet_["eSQ"])
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        ig = ig_['OBJ']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['OB'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        ig = ig_['C5']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['X1X2'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['C8']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['X1SQ'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.008))
        ig = ig_['C9']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['X2SQ'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN                -59.01312394
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.cupper[np.arange(pb.nle)] = np.zeros((pb.nle,1))
        pb.clower[np.arange(pb.nle+pb.neq,pb.m)] = np.zeros((pb.nge,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "OQR2-AN-2-5"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eOBNL(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        B3 = .1269366345
        B4 = -0.20567665
        B4 = 0.01*B4
        B5 = 0.103450e-4
        B7 = .0302344793
        B8 = -0.12813448
        B8 = 0.01*B8
        B9 = 0.352599e-4
        B10 = -0.2266e-6
        B11 = 0.2564581253
        B12 = -.003460403
        B13 = 0.135139e-4
        B14 = -.1064434908
        B14 = B14-28.0
        B15 = -0.52375e-5
        B16 = -0.63e-8
        B17 = 0.7e-9
        B18 = 0.3405462
        B18 = 0.001*B18
        B19 = -0.16638e-5
        B20 = -2.86731123
        B20 = B20-0.92e-8
        A = B7*EV_[0]+B8*EV_[0]**2+B9*EV_[0]**3+B10*EV_[0]**4
        DADX = B7+2.0*B8*EV_[0]+3.0*B9*EV_[0]**2+4.0*B10*EV_[0]**3
        D2ADXX = 2.0*B8+6.0*B9*EV_[0]+12.0*B10*EV_[0]**2
        B = B18*EV_[0]+B15*EV_[0]**2+B16*EV_[0]**3
        DBDX = B18+2.0*B15*EV_[0]+3.0*B16*EV_[0]**2
        D2BDXX = 2.0*B15+6.0*B16*EV_[0]
        C = B3*EV_[0]**2+B4*EV_[0]**3+B5*EV_[0]**4
        DCDX = 2.0*B3*EV_[0]+3.0*B4*EV_[0]**2+4.0*B5*EV_[0]**3
        D2CDXX = 2.0*B3+6.0*B4*EV_[0]+12.0*B5*EV_[0]**2
        F = B11*EV_[1]**2+B12*EV_[1]**3+B13*EV_[1]**4
        DFDY = 2.0*B11*EV_[1]+3.0*B12*EV_[1]**2+4.0*B13*EV_[1]**3
        D2FDYY = 2.0*B11+6.0*B12*EV_[1]+12.0*B13*EV_[1]**2
        G = B17*EV_[0]**3+B19*EV_[0]
        DGDX = B19+3.0*B17*EV_[0]**2
        D2GDXX = 6.0*B17*EV_[0]
        E = np.exp(0.0005*EV_[0]*EV_[1])
        DEDX = 0.0005*EV_[1]*E
        DEDY = 0.0005*EV_[0]*E
        D2EDXX = 0.0005*EV_[1]*DEDX
        D2EDXY = 0.0005*(EV_[1]*DEDY+E)
        D2EDYY = 0.0005*EV_[0]*DEDY
        f_   = C+EV_[1]*A+F+B14/(1.0+EV_[1])+B*EV_[1]**2+G*EV_[1]**3+B20*E
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = DCDX+EV_[1]*DADX+DBDX*EV_[1]**2+DGDX*EV_[1]**3+B20*DEDX
            g_[1] = A+DFDY-B14/(1.0+EV_[1])**2+2.0*B*EV_[1]+3.0*G*EV_[1]**2+B20*DEDY
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = D2CDXX+EV_[1]*D2ADXX+D2BDXX*EV_[1]**2+D2GDXX*EV_[1]**3+B20*D2EDXX
                H_[0,1] = DADX+2.0*EV_[1]*DBDX+3.0*DGDX*EV_[1]**2+B20*D2EDXY
                H_[1,0] = H_[0,1]
                H_[1,1] = D2FDYY+2.0*B14/(1.0+EV_[1])**3+2.0*B+6.0*G*EV_[1]+B20*D2EDYY
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def en2PR(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]
            g_[1] = EV_[0]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = 1.0
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eSQ(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[0]+EV_[0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

