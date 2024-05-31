from s2xlib import *
class  CHENHARK(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CHENHARK
#    --------
# 
#    A bound-constrained version the Linear Complementarity problem
# 
#    Find x such that w = M x + q, x and w nonnegative and x^T w = 0,
#    where
# 
#    M = (  6   -4   1   0  ........ 0 ) 
#        ( -4    6  -4   1  ........ 0 )
#        (  1   -4   6  -4  ........ 0 )
#        (  0    1  -4   6  ........ 0 )  
#           ..........................
#        (  0   ........... 0  1 -4  6 )
# 
#    and q is given.
# 
#    Source: 
#    B. Chen and P. T. Harker,
#    SIMAX 14 (1993) 1168-1190
# 
#    SDIF input: Nick Gould, November 1993.
# 
#    classification = "QBR2-AN-V-V"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'CHENHARK'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'CHENHARK'
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
# IE N                   100            $-PARAMETER
# IE N                   1000           $-PARAMETER     original value
# IE N                   5000           $-PARAMETER
# IE N                   10000          $-PARAMETER
# IE N                   50000          $-PARAMETER
        if nargin<2:
            v_['NFREE'] = int(5);  #  SIF file default value
        else:
            v_['NFREE'] = int(args[1])
# IE NFREE               50             $-PARAMETER
# IE NFREE               500            $-PARAMETER     original value
# IE NFREE               2500           $-PARAMETER
# IE NFREE               5000           $-PARAMETER
# IE NFREE               10000          $-PARAMETER
        if nargin<3:
            v_['NDEGEN'] = int(2);  #  SIF file default value
        else:
            v_['NDEGEN'] = int(args[2])
# IE NDEGEN              20             $-PARAMETER
# IE NDEGEN              200            $-PARAMETER     original value
# IE NDEGEN              500            $-PARAMETER
# IE NDEGEN              1000           $-PARAMETER
# IE NDEGEN              2000           $-PARAMETER
        v_['-1'] = -1
        v_['0'] = 0
        v_['1'] = 1
        v_['2'] = 2
        v_['N-1'] = -1+v_['N']
        v_['N+1'] = 1+v_['N']
        v_['N+2'] = 2+v_['N']
        v_['NFREE+1'] = 1+v_['NFREE']
        v_['NF+ND'] = v_['NFREE']+v_['NDEGEN']
        v_['NF+ND+1'] = 1+v_['NF+ND']
        v_['X'+str(int(v_['-1']))] = 0.0
        v_['X'+str(int(v_['0']))] = 0.0
        for I in range(int(v_['1']),int(v_['NFREE'])+1):
            v_['X'+str(I)] = 1.0
        for I in range(int(v_['NFREE+1']),int(v_['N+2'])+1):
            v_['X'+str(I)] = 0.0
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
        for I in range(int(v_['2']),int(v_['N-1'])+1):
            v_['I+1'] = 1+I
            v_['I-1'] = -1+I
            [ig,ig_,_] = s2x_ii('Q'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(int(v_['I+1']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['I-1']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = float(-2.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('Q'+str(int(v_['0'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X'+str(int(v_['1']))]
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('Q'+str(int(v_['1'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X'+str(int(v_['1']))]
        pbm.A[ig,iv] = float(2.0)+pbm.A[ig,iv]
        iv = ix_['X'+str(int(v_['2']))]
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('Q'+str(int(v_['N'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X'+str(int(v_['N']))]
        pbm.A[ig,iv] = float(2.0)+pbm.A[ig,iv]
        iv = ix_['X'+str(int(v_['N-1']))]
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('Q'+str(int(v_['N+1'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X'+str(int(v_['N']))]
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['NF+ND'])+1):
            v_['I+1'] = 1+I
            v_['I+2'] = 2+I
            v_['I-1'] = -1+I
            v_['I-2'] = -2+I
            v_['Q1'] = -6.0*v_['X'+str(I)]
            v_['Q2'] = 4.0*v_['X'+str(int(v_['I+1']))]
            v_['Q3'] = 4.0*v_['X'+str(int(v_['I-1']))]
            v_['Q4'] = -1.0*v_['X'+str(int(v_['I+2']))]
            v_['Q5'] = -1.0*v_['X'+str(int(v_['I-2']))]
            v_['Q'] = v_['Q1']+v_['Q2']
            v_['Q'] = v_['Q']+v_['Q3']
            v_['Q'] = v_['Q']+v_['Q4']
            v_['Q'] = v_['Q']+v_['Q5']
            [ig,ig_,_] = s2x_ii('L',ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = float(v_['Q'])+pbm.A[ig,iv]
        for I in range(int(v_['NF+ND+1']),int(v_['N'])+1):
            v_['I+1'] = 1+I
            v_['I+2'] = 2+I
            v_['I-1'] = -1+I
            v_['I-2'] = -2+I
            v_['Q1'] = -6.0*v_['X'+str(I)]
            v_['Q2'] = 4.0*v_['X'+str(int(v_['I+1']))]
            v_['Q3'] = 4.0*v_['X'+str(int(v_['I-1']))]
            v_['Q4'] = -1.0*v_['X'+str(int(v_['I+2']))]
            v_['Q5'] = -1.0*v_['X'+str(int(v_['I-2']))]
            v_['Q'] = v_['Q1']+v_['Q2']
            v_['Q'] = v_['Q']+v_['Q3']
            v_['Q'] = v_['Q']+v_['Q4']
            v_['Q'] = v_['Q']+v_['Q5']
            v_['Q'] = 1.0+v_['Q']
            [ig,ig_,_] = s2x_ii('L',ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = float(v_['Q'])+pbm.A[ig,iv]
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        for I in range(int(v_['1']),int(v_['N'])+1):
            pb.x0[ix_['X'+str(I)]] = float(0.5)
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
        for I in range(int(v_['0']),int(v_['N+1'])+1):
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

