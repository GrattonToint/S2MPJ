from s2xlib import *
class  HARKERP2(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HARKERP2
#    --------
# 
#    A bound-constrained version of a Linear Complementarity problem
#    posed by Harker and Pang.
# 
#    Source: 
#    P. T. Harker and J.-S. Pang,
#    "A damped Newton method for the linear complementarity problem",
#    in 'Allgower and Georg: Computational solution of nonlinear
#    systems of equations', AMS lectures in Applied Mathematics 26,
#    AMS, Providence, Rhode Island, USA, pp 265-284.
# 
#    SIF input: Nick Gould, July 1993.
# 
#    classification = "QBR2-AN-V-V"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HARKERP2'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'HARKERP2'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(10);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
#           Alternative values for the SIF file parameters:
# IE N                   100            $-PARAMETER     original value
# IE N                   500            $-PARAMETER
# IE N                   1000           $-PARAMETER
# IE N                   5000           $-PARAMETER
# IE N                   10000          $-PARAMETER
        v_['0'] = 0
        v_['1'] = 1
        v_['2'] = 2
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
        for I in range(int(v_['1']),int(v_['N'])+1):
            [ig,ig_,_] = s2x_ii('S'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            pbm.gscale = arrset(pbm.gscale,ig,float(-1.0))
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('Q'+str(int(v_['0'])),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('Q'+str(int(v_['1'])),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('Q'+str(int(v_['1'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        pbm.gscale = arrset(pbm.gscale,ig,float(0.5))
        for J in range(int(v_['2']),int(v_['N'])+1):
            [ig,ig_,_] = s2x_ii('Q'+str(J),ig_)
            gtype = arrset(gtype,ig,'<>')
            pbm.gscale = arrset(pbm.gscale,ig,float(0.25))
            for I in range(int(J),int(v_['N'])+1):
                [ig,ig_,_] = s2x_ii('Q'+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
                iv = ix_['X'+str(I)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['RI'] = float(I)
            pb.x0[ix_['X'+str(I)]] = float(v_['RI'])
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2x_ii('gHALFL2',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            ig = ig_['S'+str(I)]
            pbm.grftype = arrset(pbm.grftype,ig,'gHALFL2')
            ig = ig_['Q'+str(I)]
            pbm.grftype = arrset(pbm.grftype,ig,'gHALFL2')
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 1.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "QBR2-AN-V-V"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gHALFL2(pbm,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= 5.0e-1*GVAR_*GVAR_
        if nargout>1:
            g_ = GVAR_
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 1.0e+0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

