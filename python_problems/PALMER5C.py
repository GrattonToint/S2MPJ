from s2mpjlib import *
class  PALMER5C(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : PALMER5C
#    *********
# 
#    A linear least squares problem
#    arising from chemical kinetics.
# 
#     model: H-N=C=Se TZVP + MP2
#    fitting Y to A0 T_0 + A2 T_2 + A4 T_4 + A6 T_6 + A8 T_8 +
#                 A10 T_10 + A12 T_12 + A14 T_14
#    where T_i is the i-th (shifted) Chebyshev polynomial
# 
#    Source:
#    M. Palmer, Edinburgh, private communication.
# 
#    SIF input: Nick Gould, 1992.
# 
#    classification = "QUR2-RN-6-0"
# 
#    Number of data points
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'PALMER5C'

    def __init__(self, *args): 
        import numpy as np
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['M'] = 23
        v_['1'] = 1
        v_['2'] = 2
        v_['12'] = 12
        v_['14'] = 14
        v_['X12'] = 0.000000
        v_['X13'] = 1.570796
        v_['X14'] = 1.396263
        v_['X15'] = 1.308997
        v_['X16'] = 1.221730
        v_['X17'] = 1.125835
        v_['X18'] = 1.047198
        v_['X19'] = 0.872665
        v_['X20'] = 0.698132
        v_['X21'] = 0.523599
        v_['X22'] = 0.349066
        v_['X23'] = 0.174533
        v_['B'] = v_['X13']
        v_['A'] = -1.0e+0*v_['B']
        v_['DIFF'] = 2.0e+0*v_['B']
        v_['Y12'] = 83.57418
        v_['Y13'] = 81.007654
        v_['Y14'] = 18.983286
        v_['Y15'] = 8.051067
        v_['Y16'] = 2.044762
        v_['Y17'] = 0.000000
        v_['Y18'] = 1.170451
        v_['Y19'] = 10.479881
        v_['Y20'] = 25.785001
        v_['Y21'] = 44.126844
        v_['Y22'] = 62.822177
        v_['Y23'] = 77.719674
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2mpj_ii('A0',ix_)
        self.xnames=arrset(self.xnames,iv,'A0')
        [iv,ix_,_] = s2mpj_ii('A2',ix_)
        self.xnames=arrset(self.xnames,iv,'A2')
        [iv,ix_,_] = s2mpj_ii('A4',ix_)
        self.xnames=arrset(self.xnames,iv,'A4')
        [iv,ix_,_] = s2mpj_ii('A6',ix_)
        self.xnames=arrset(self.xnames,iv,'A6')
        [iv,ix_,_] = s2mpj_ii('A8',ix_)
        self.xnames=arrset(self.xnames,iv,'A8')
        [iv,ix_,_] = s2mpj_ii('A10',ix_)
        self.xnames=arrset(self.xnames,iv,'A10')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.A       = lil_matrix((1000000,1000000))
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames      = np.array([])
        self.cnames = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['12']),int(v_['M'])+1):
            v_['T0'] = 1.0e+0
            v_['Y'] = 2.0e+0*v_['X'+str(I)]
            v_['Y'] = v_['Y']-v_['A']
            v_['Y'] = v_['Y']-v_['B']
            v_['Y'] = v_['Y']/v_['DIFF']
            v_['T1'] = v_['Y']
            v_['2Y'] = 2.0e+0*v_['Y']
            for J in range(int(v_['2']),int(v_['14'])+1):
                v_['J-1'] = -1+J
                v_['J-2'] = -2+J
                v_['T'+str(J)] = v_['2Y']*v_['T'+str(int(v_['J-1']))]
                v_['T'+str(J)] = v_['T'+str(J)]-v_['T'+str(int(v_['J-2']))]
            [ig,ig_,_] = s2mpj_ii('O'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['A0']
            self.A[ig,iv] = float(v_['T0'])+self.A[ig,iv]
            iv = ix_['A2']
            self.A[ig,iv] = float(v_['T2'])+self.A[ig,iv]
            iv = ix_['A4']
            self.A[ig,iv] = float(v_['T4'])+self.A[ig,iv]
            iv = ix_['A6']
            self.A[ig,iv] = float(v_['T6'])+self.A[ig,iv]
            iv = ix_['A8']
            self.A[ig,iv] = float(v_['T8'])+self.A[ig,iv]
            iv = ix_['A10']
            self.A[ig,iv] = float(v_['T10'])+self.A[ig,iv]
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        for I in range(int(v_['12']),int(v_['M'])+1):
            self.gconst = arrset(self.gconst,ig_['O'+str(I)],float(v_['Y'+str(I)]))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        self.xlower[ix_['A0']] = -float('Inf')
        self.xupper[ix_['A0']] = +float('Inf')
        self.xlower[ix_['A2']] = -float('Inf')
        self.xupper[ix_['A2']] = +float('Inf')
        self.xlower[ix_['A4']] = -float('Inf')
        self.xupper[ix_['A4']] = +float('Inf')
        self.xlower[ix_['A6']] = -float('Inf')
        self.xupper[ix_['A6']] = +float('Inf')
        self.xlower[ix_['A8']] = -float('Inf')
        self.xupper[ix_['A8']] = +float('Inf')
        self.xlower[ix_['A10']] = -float('Inf')
        self.xupper[ix_['A10']] = +float('Inf')
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(1.0))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gL2',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['12']),int(v_['M'])+1):
            ig = ig_['O'+str(I)]
            self.grftype = arrset(self.grftype,ig,'gL2')
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        self.objlower = 0.0
#    Solution
# LO SOLTN              5.0310687D-02
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        self.A.resize(ngrp,self.n)
        self.A     = self.A.tocsr()
        sA1,sA2    = self.A.shape
        self.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass = "QUR2-RN-6-0"
# ********************
#  SET UP THE GROUPS *
#  ROUTINE           *
# ********************

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gL2(self,nargout,*args):

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

