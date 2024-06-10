from s2mpjlib import *
class  CURLY30(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CURLY30
#    --------
# 
#    A banded function with semi-bandwidth 30 and
#    negative curvature near the starting point
# 
#    Source: Nick Gould
# 
#    SIF input: Nick Gould, September 1997.
# 
#    classification = "OUR2-AN-V-0"
# 
#    Number of variables
# 
#           Alternative values for the SIF file parameters:
# IE N                   100            $-PARAMETER
# IE N                   1000           $-PARAMETER     original value
# IE N                   10000          $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'CURLY30'

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
        if nargin<1:
            v_['N'] = int(35);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
        v_['K'] = 30
        v_['1'] = 1
        v_['N-K'] = v_['N']-v_['K']
        v_['N-K+1'] = 1+v_['N-K']
        v_['RN'] = float(v_['N'])
        v_['RN+1'] = 1+v_['RN']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['N-K'])+1):
            v_['I+K'] = I+v_['K']
            for J in range(int(I),int(v_['I+K'])+1):
                [ig,ig_,_] = s2mpj_ii('Q'+str(I),ig_)
                gtype = arrset(gtype,ig,'<>')
                iv = ix_['X'+str(J)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        for I in range(int(v_['N-K+1']),int(v_['N'])+1):
            for J in range(int(I),int(v_['N'])+1):
                [ig,ig_,_] = s2mpj_ii('Q'+str(I),ig_)
                gtype = arrset(gtype,ig,'<>')
                iv = ix_['X'+str(J)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower = np.zeros((pb.n,1))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['RI'] = float(I)
            v_['T'] = v_['RI']/v_['RN+1']
            v_['T'] = 0.0001*v_['T']
            pb.x0[ix_['X'+str(I)]] = float(v_['T'])
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gP4',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for ig in range(0,ngrp):
            pbm.grftype = arrset(pbm.grftype,ig,'gP4')
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# ZL SOLTN               -1.003163D+5   $ (n=1000)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "OUR2-AN-V-0"
        self.pb = pb; self.pbm = pbm
# ********************
#  SET UP THE GROUPS *
#  ROUTINE           *
# ********************

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gP4(pbm,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        APB = 2.0e+1
        f_= GVAR_*(GVAR_*(GVAR_**2-APB)-1.0e-1)
        if nargout>1:
            g_ = 2.0e+0*GVAR_*(2.0e+0*GVAR_**2-APB)-1.0e-1
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 1.2e+1*GVAR_**2-2.0e+0*APB
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

