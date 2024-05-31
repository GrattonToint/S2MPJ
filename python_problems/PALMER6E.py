from s2xlib import *
class  PALMER6E(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : PALMER6E
#    *********
# 
#    A nonlinear least squares problem
#    arising from chemical kinetics.
# 
#    model: H-N=C=Se TZVP + MP2
#    fitting Y to A0 + A2 X**2 + A4 X**4 + A6 X**6 + A8 X**8 +
#                 A10 X**10 + L * EXP( -K X**2 )
# 
#    Source:
#    M.  Palmer, Edinburgh, private communication.
# 
#    SIF input: Nick Gould, 1992.
# 
#    classification = "SBR2-RN-8-0"
# 
#    Number of data points
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'PALMER6E'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'PALMER6E'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['M'] = 24
        v_['1'] = 1
        v_['12'] = 12
        v_['X12'] = 0.000000
        v_['X13'] = 1.570796
        v_['X14'] = 1.396263
        v_['X15'] = 1.221730
        v_['X16'] = 1.047198
        v_['X17'] = 0.872665
        v_['X18'] = 0.785398
        v_['X19'] = 0.732789
        v_['X20'] = 0.698132
        v_['X21'] = 0.610865
        v_['X22'] = 0.523599
        v_['X23'] = 0.349066
        v_['X24'] = 0.174533
        v_['Y12'] = 10.678659
        v_['Y13'] = 75.414511
        v_['Y14'] = 41.513459
        v_['Y15'] = 20.104735
        v_['Y16'] = 7.432436
        v_['Y17'] = 1.298082
        v_['Y18'] = 0.171300
        v_['Y19'] = 0.000000
        v_['Y20'] = 0.068203
        v_['Y21'] = 0.774499
        v_['Y22'] = 2.070002
        v_['Y23'] = 5.574556
        v_['Y24'] = 9.026378
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2x_ii('A0',ix_)
        pb.xnames=arrset(pb.xnames,iv,'A0')
        [iv,ix_,_] = s2x_ii('A2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'A2')
        [iv,ix_,_] = s2x_ii('A4',ix_)
        pb.xnames=arrset(pb.xnames,iv,'A4')
        [iv,ix_,_] = s2x_ii('A6',ix_)
        pb.xnames=arrset(pb.xnames,iv,'A6')
        [iv,ix_,_] = s2x_ii('A8',ix_)
        pb.xnames=arrset(pb.xnames,iv,'A8')
        [iv,ix_,_] = s2x_ii('A10',ix_)
        pb.xnames=arrset(pb.xnames,iv,'A10')
        [iv,ix_,_] = s2x_ii('K',ix_)
        pb.xnames=arrset(pb.xnames,iv,'K')
        [iv,ix_,_] = s2x_ii('L',ix_)
        pb.xnames=arrset(pb.xnames,iv,'L')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['12']),int(v_['M'])+1):
            v_['XSQR'] = v_['X'+str(I)]*v_['X'+str(I)]
            v_['XQUART'] = v_['XSQR']*v_['XSQR']
            v_['X**6'] = v_['XSQR']*v_['XQUART']
            v_['X**8'] = v_['XSQR']*v_['X**6']
            v_['X**10'] = v_['XSQR']*v_['X**8']
            v_['X**12'] = v_['XSQR']*v_['X**10']
            v_['X**14'] = v_['XSQR']*v_['X**12']
            [ig,ig_,_] = s2x_ii('O'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['A0']
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['A2']
            pbm.A[ig,iv] = float(v_['XSQR'])+pbm.A[ig,iv]
            iv = ix_['A4']
            pbm.A[ig,iv] = float(v_['XQUART'])+pbm.A[ig,iv]
            iv = ix_['A6']
            pbm.A[ig,iv] = float(v_['X**6'])+pbm.A[ig,iv]
            iv = ix_['A8']
            pbm.A[ig,iv] = float(v_['X**8'])+pbm.A[ig,iv]
            iv = ix_['A10']
            pbm.A[ig,iv] = float(v_['X**10'])+pbm.A[ig,iv]
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        for I in range(int(v_['12']),int(v_['M'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['O'+str(I)],float(v_['Y'+str(I)]))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('inf'))
        pb.xupper = np.full((pb.n,1),float('inf'))
        pb.xlower[ix_['A0']] = -float('Inf')
        pb.xupper[ix_['A0']] = +float('Inf')
        pb.xlower[ix_['A2']] = -float('Inf')
        pb.xupper[ix_['A2']] = +float('Inf')
        pb.xlower[ix_['A4']] = -float('Inf')
        pb.xupper[ix_['A4']] = +float('Inf')
        pb.xlower[ix_['A6']] = -float('Inf')
        pb.xupper[ix_['A6']] = +float('Inf')
        pb.xlower[ix_['A8']] = -float('Inf')
        pb.xupper[ix_['A8']] = +float('Inf')
        pb.xlower[ix_['A10']] = -float('Inf')
        pb.xupper[ix_['A10']] = +float('Inf')
        pb.xlower[ix_['L']] = -float('Inf')
        pb.xupper[ix_['L']] = +float('Inf')
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(1.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'ePROD', iet_)
        elftv = loaset(elftv,it,0,'K')
        elftv = loaset(elftv,it,1,'L')
        elftp = []
        elftp = loaset(elftp,it,0,'XSQR')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['12']),int(v_['M'])+1):
            v_['XSQR'] = v_['X'+str(I)]*v_['X'+str(I)]
            ename = 'E'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'ePROD')
            ielftype = arrset(ielftype, ie, iet_["ePROD"])
            vname = 'K'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,1.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='K')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'L'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,1.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='L')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='XSQR')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['XSQR']))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2x_ii('gL2',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['12']),int(v_['M'])+1):
            ig = ig_['O'+str(I)]
            pbm.grftype = arrset(pbm.grftype,ig,'gL2')
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "SBR2-RN-8-0"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def ePROD(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        EXPON = np.exp(-EV_[0]*pbm.elpar[iel_][0])
        f_   = EV_[1]*EXPON
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -pbm.elpar[iel_][0]*EV_[1]*EXPON
            g_[1] = EXPON
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = pbm.elpar[iel_][0]*pbm.elpar[iel_][0]*EV_[1]*EXPON
                H_[0,1] = -pbm.elpar[iel_][0]*EXPON
                H_[1,0] = H_[0,1]
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

