from s2mpjlib import *
class  STNQP1(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : STNQP1
#    *********
# 
#    A non-convex quadratic program with some structure.
# 
#    The objective function is of the form
#       sum (i=0,n) x_i^2 - 0.5 sum (l=1,n/p) sum(i=1,p) sum(k;i) x_{k+l}^2,
#    where n = 2^p and (k;i) means k takes the values of the first i powers of 2
#    eg, (k:3) = {k = {1,2,4}} and (k:7) = {k = {1,2,4,8,16,32}}.
#    There are equality constraints of the form
#    
#       sum(k;i) x_{k+l-1} = i, where l=1,n/p,2 and i=1,p.
#    Finally, there are simple bounds
#          2 <= x_i, y_i <= 2    (i=0,n).
# 
#    SIF input: Nick Gould, May 1996
# 
#    classification = "QLR2-AN-V-V"
# 
#    There will be 2**p + 1 variables
# 
#           Alternative values for the SIF file parameters:
# IE P                   2              $-PARAMETER n = 5
# IE P                   4              $-PARAMETER n = 17
# IE P                   6              $-PARAMETER n = 65
# IE P                   8              $-PARAMETER n = 257
# IE P                   10             $-PARAMETER n = 1025
# IE P                   12             $-PARAMETER n = 4097     original value
# IE P                   13             $-PARAMETER n = 8193
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'STNQP1'

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
            v_['P'] = int(4);  #  SIF file default value
        else:
            v_['P'] = int(args[0])
# IE P                   14             $-PARAMETER n = 16395
# IE P                   15             $-PARAMETER n = 32769
# IE P                   16             $-PARAMETER n = 65537
        v_['0'] = 0
        v_['1'] = 1
        v_['2'] = 2
        v_['N'] = 1
        for I in range(int(v_['1']),int(v_['P'])+1):
            v_['N'] = v_['N']*v_['2']
        v_['N/P'] = int(np.fix(v_['N']/v_['P']))
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['0']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['0']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('O'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        for L in range(int(v_['1']),int(v_['N/P'])+1):
            for I in range(int(v_['1']),int(v_['P'])+1):
                v_['K'] = v_['1']
                for J in range(int(v_['1']),int(I)+1):
                    v_['K+L'] = v_['K']+L
                    [ig,ig_,_] = s2mpj_ii('N'+str(I)+','+str(L),ig_)
                    gtype = arrset(gtype,ig,'<>')
                    iv = ix_['X'+str(int(v_['K+L']))]
                    pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
                    v_['K'] = v_['K']*v_['2']
        for L in range(int(v_['1']),int(v_['N/P'])+1,int(v_['2'])):
            for I in range(int(v_['1']),int(v_['P'])+1):
                v_['K'] = v_['1']
                for J in range(int(v_['1']),int(I)+1):
                    v_['K-1'] = v_['K']-v_['1']
                    v_['K+L-1'] = v_['K-1']+L
                    [ig,ig_,_] = s2mpj_ii('E'+str(I)+','+str(L),ig_)
                    gtype = arrset(gtype,ig,'==')
                    cnames = arrset(cnames,ig,'E'+str(I)+','+str(L))
                    iv = ix_['X'+str(int(v_['K+L-1']))]
                    pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
                    v_['K'] = v_['K']*v_['2']
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
        for L in range(int(v_['1']),int(v_['N/P'])+1,int(v_['2'])):
            for I in range(int(v_['1']),int(v_['P'])+1):
                v_['RI'] = float(I)
                pbm.gconst = arrset(pbm.gconst,ig_['E'+str(I)+','+str(L)],float(v_['RI']))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),float('inf'))
        for I in range(int(v_['0']),int(v_['N'])+1):
            pb.xlower[ix_['X'+str(I)]] = -2.0
            pb.xupper[ix_['X'+str(I)]] = 2.0
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(0.5))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gPSQR',igt_)
        [it,igt_,_] = s2mpj_ii('gPSQR',igt_)
        grftp = []
        grftp = loaset(grftp,it,0,'P')
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        pbm.grpar   = []
        for I in range(int(v_['0']),int(v_['N'])+1):
            ig = ig_['O'+str(I)]
            pbm.grftype = arrset(pbm.grftype,ig,'gPSQR')
            posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P')
            pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(1.0))
        for L in range(int(v_['1']),int(v_['N/P'])+1):
            for I in range(int(v_['1']),int(v_['P'])+1):
                ig = ig_['N'+str(I)+','+str(L)]
                pbm.grftype = arrset(pbm.grftype,ig,'gPSQR')
                posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P')
                pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(-0.5))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLUTION            -1.361565E+5   $ (P=12)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.lincons   = np.arange(len(pbm.congrps))
        pb.pbclass = "QLR2-AN-V-V"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gPSQR(pbm,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= pbm.grpar[igr_][0]*GVAR_*GVAR_
        if nargout>1:
            g_ = 2.0*pbm.grpar[igr_][0]*GVAR_
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 2.0*pbm.grpar[igr_][0]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

