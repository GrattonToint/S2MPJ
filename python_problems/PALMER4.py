from s2xlib import *
class  PALMER4(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : PALMER4
#    *********
# 
#    A nonlinear least squares problem with bounds
#    arising from chemical kinetics.
# 
#    model: H-N=C=Se TZVP + MP2
#    fitting Y to A X**2 + B / ( C + X**2 / D ), B, C, D nonnegative.
# 
#    Source:
#    M. Palmer, Edinburgh, private communication.
# 
#    SIF input: Nick Gould, 1990.
# 
#    classification = "SBR2-RN-4-0"
# 
#    Number of data points
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'PALMER4'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'PALMER4'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['M'] = 23
        v_['1'] = 1
        v_['X1'] = -1.658063
        v_['X2'] = -1.570796
        v_['X3'] = -1.396263
        v_['X4'] = -1.221730
        v_['X5'] = -1.047198
        v_['X6'] = -0.872665
        v_['X7'] = -0.741119
        v_['X8'] = -0.698132
        v_['X9'] = -0.523599
        v_['X10'] = -0.349066
        v_['X11'] = -0.174533
        v_['X12'] = 0.0
        v_['X13'] = 0.174533
        v_['X14'] = 0.349066
        v_['X15'] = 0.523599
        v_['X16'] = 0.698132
        v_['X17'] = 0.741119
        v_['X18'] = 0.872665
        v_['X19'] = 1.047198
        v_['X20'] = 1.221730
        v_['X21'] = 1.396263
        v_['X22'] = 1.570796
        v_['X23'] = 1.658063
        v_['Y1'] = 67.27625
        v_['Y2'] = 52.8537
        v_['Y3'] = 30.2718
        v_['Y4'] = 14.9888
        v_['Y5'] = 5.5675
        v_['Y6'] = 0.92603
        v_['Y7'] = 0.0
        v_['Y8'] = 0.085108
        v_['Y9'] = 1.867422
        v_['Y10'] = 5.014768
        v_['Y11'] = 8.263520
        v_['Y12'] = 9.8046208
        v_['Y13'] = 8.263520
        v_['Y14'] = 5.014768
        v_['Y15'] = 1.867422
        v_['Y16'] = 0.085108
        v_['Y17'] = 0.0
        v_['Y18'] = 0.92603
        v_['Y19'] = 5.5675
        v_['Y20'] = 14.9888
        v_['Y21'] = 30.2718
        v_['Y22'] = 52.8537
        v_['Y23'] = 67.27625
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2x_ii('A',ix_)
        pb.xnames=arrset(pb.xnames,iv,'A')
        [iv,ix_,_] = s2x_ii('B',ix_)
        pb.xnames=arrset(pb.xnames,iv,'B')
        [iv,ix_,_] = s2x_ii('C',ix_)
        pb.xnames=arrset(pb.xnames,iv,'C')
        [iv,ix_,_] = s2x_ii('D',ix_)
        pb.xnames=arrset(pb.xnames,iv,'D')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['M'])+1):
            v_['XSQR'] = v_['X'+str(I)]*v_['X'+str(I)]
            [ig,ig_,_] = s2x_ii('O'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['A']
            pbm.A[ig,iv] = float(v_['XSQR'])+pbm.A[ig,iv]
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        for I in range(int(v_['1']),int(v_['M'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['O'+str(I)],float(v_['Y'+str(I)]))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('inf'))
        pb.xupper = np.full((pb.n,1),float('inf'))
        pb.xlower[ix_['A']] = -float('Inf')
        pb.xupper[ix_['A']] = +float('Inf')
        pb.xlower[ix_['B']] = 0.00001
        pb.xlower[ix_['C']] = 0.00001
        pb.xlower[ix_['D']] = 0.00001
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(1.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eQUOT', iet_)
        elftv = loaset(elftv,it,0,'B')
        elftv = loaset(elftv,it,1,'C')
        elftv = loaset(elftv,it,2,'D')
        elftp = []
        elftp = loaset(elftp,it,0,'XSQR')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['1']),int(v_['M'])+1):
            v_['XSQR'] = v_['X'+str(I)]*v_['X'+str(I)]
            ename = 'E'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eQUOT')
            ielftype = arrset(ielftype, ie, iet_["eQUOT"])
            vname = 'B'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,1.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='B')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'C'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,1.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='C')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'D'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,1.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='D')
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
        for I in range(int(v_['1']),int(v_['M'])+1):
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
        pb.pbclass = "SBR2-RN-4-0"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eQUOT(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        DENOM = 1.0/(EV_[1]+pbm.elpar[iel_][0]/EV_[2])
        f_   = EV_[0]*DENOM
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = DENOM
            g_[1] = -EV_[0]*DENOM*DENOM
            g_[2] = EV_[0]*pbm.elpar[iel_][0]*(DENOM/EV_[2])**2
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = -DENOM*DENOM
                H_[1,0] = H_[0,1]
                H_[0,2] = pbm.elpar[iel_][0]*(DENOM/EV_[2])**2
                H_[2,0] = H_[0,2]
                H_[1,1] = 2.0*EV_[0]*DENOM**3
                H_[1,2] = -2.0*EV_[0]*pbm.elpar[iel_][0]*DENOM*(DENOM/EV_[2])**2
                H_[2,1] = H_[1,2]
                H_[2,2]  = (
                      2.0*EV_[0]*DENOM**3*(pbm.elpar[iel_][0]/EV_[2]**2)**2-2.0*EV_[0]*DENOM**2*pbm.elpar[iel_][0]/EV_[2]**3)
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
