from s2mpjlib import *
class  PALMER1D(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : PALMER1D
#    *********
# 
#    A linear least squares problem
#    arising from chemical kinetics.
# 
#    model: H-N=N=N TZVP+MP2
#    fitting Y to A0 + A2 X**2 + A4 X**4 + A6 X**6 + A8 X**8 +
#                 A10 X**10 + A12 X**12
# 
#    Source:
#    M. Palmer, Edinburgh, private communication.
# 
#    SIF input: Nick Gould, 1990.
# 
#    classification = "QUR2-RN-7-0"
# 
#    Number of data points
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'PALMER1D'

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
        v_['M'] = 35
        v_['1'] = 1
        v_['X1'] = -1.788963
        v_['X2'] = -1.745329
        v_['X3'] = -1.658063
        v_['X4'] = -1.570796
        v_['X5'] = -1.483530
        v_['X6'] = -1.396263
        v_['X7'] = -1.308997
        v_['X8'] = -1.218612
        v_['X9'] = -1.134464
        v_['X10'] = -1.047198
        v_['X11'] = -0.872665
        v_['X12'] = -0.698132
        v_['X13'] = -0.523599
        v_['X14'] = -0.349066
        v_['X15'] = -0.174533
        v_['X16'] = 0.0000000
        v_['X17'] = 1.788963
        v_['X18'] = 1.745329
        v_['X19'] = 1.658063
        v_['X20'] = 1.570796
        v_['X21'] = 1.483530
        v_['X22'] = 1.396263
        v_['X23'] = 1.308997
        v_['X24'] = 1.218612
        v_['X25'] = 1.134464
        v_['X26'] = 1.047198
        v_['X27'] = 0.872665
        v_['X28'] = 0.698132
        v_['X29'] = 0.523599
        v_['X30'] = 0.349066
        v_['X31'] = 0.174533
        v_['X32'] = -1.8762289
        v_['X33'] = -1.8325957
        v_['X34'] = 1.8762289
        v_['X35'] = 1.8325957
        v_['Y1'] = 78.596218
        v_['Y2'] = 65.77963
        v_['Y3'] = 43.96947
        v_['Y4'] = 27.038816
        v_['Y5'] = 14.6126
        v_['Y6'] = 6.2614
        v_['Y7'] = 1.538330
        v_['Y8'] = 0.000000
        v_['Y9'] = 1.188045
        v_['Y10'] = 4.6841
        v_['Y11'] = 16.9321
        v_['Y12'] = 33.6988
        v_['Y13'] = 52.3664
        v_['Y14'] = 70.1630
        v_['Y15'] = 83.4221
        v_['Y16'] = 88.3995
        v_['Y17'] = 78.596218
        v_['Y18'] = 65.77963
        v_['Y19'] = 43.96947
        v_['Y20'] = 27.038816
        v_['Y21'] = 14.6126
        v_['Y22'] = 6.2614
        v_['Y23'] = 1.538330
        v_['Y24'] = 0.000000
        v_['Y25'] = 1.188045
        v_['Y26'] = 4.6841
        v_['Y27'] = 16.9321
        v_['Y28'] = 33.6988
        v_['Y29'] = 52.3664
        v_['Y30'] = 70.1630
        v_['Y31'] = 83.4221
        v_['Y32'] = 108.18086
        v_['Y33'] = 92.733676
        v_['Y34'] = 108.18086
        v_['Y35'] = 92.733676
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2mpj_ii('A0',ix_)
        pb.xnames=arrset(pb.xnames,iv,'A0')
        [iv,ix_,_] = s2mpj_ii('A2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'A2')
        [iv,ix_,_] = s2mpj_ii('A4',ix_)
        pb.xnames=arrset(pb.xnames,iv,'A4')
        [iv,ix_,_] = s2mpj_ii('A6',ix_)
        pb.xnames=arrset(pb.xnames,iv,'A6')
        [iv,ix_,_] = s2mpj_ii('A8',ix_)
        pb.xnames=arrset(pb.xnames,iv,'A8')
        [iv,ix_,_] = s2mpj_ii('A10',ix_)
        pb.xnames=arrset(pb.xnames,iv,'A10')
        [iv,ix_,_] = s2mpj_ii('A12',ix_)
        pb.xnames=arrset(pb.xnames,iv,'A12')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['M'])+1):
            v_['XSQR'] = v_['X'+str(I)]*v_['X'+str(I)]
            v_['XQUART'] = v_['XSQR']*v_['XSQR']
            v_['X**6'] = v_['XSQR']*v_['XQUART']
            v_['X**8'] = v_['XSQR']*v_['X**6']
            v_['X**10'] = v_['XSQR']*v_['X**8']
            v_['X**12'] = v_['XSQR']*v_['X**10']
            v_['X**14'] = v_['XSQR']*v_['X**12']
            [ig,ig_,_] = s2mpj_ii('O'+str(I),ig_)
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
            iv = ix_['A12']
            pbm.A[ig,iv] = float(v_['X**12'])+pbm.A[ig,iv]
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        for I in range(int(v_['1']),int(v_['M'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['O'+str(I)],float(v_['Y'+str(I)]))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
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
        pb.xlower[ix_['A12']] = -float('Inf')
        pb.xupper[ix_['A12']] = +float('Inf')
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(1.0))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gL2',igt_)
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
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        pb.objlower = 0.0
#    Solution
# LO SOLTN              0.652673985
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "QUR2-RN-7-0"
        self.pb = pb; self.pbm = pbm
# ********************
#  SET UP THE GROUPS *
#  ROUTINE           *
# ********************

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

