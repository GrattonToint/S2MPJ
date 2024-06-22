from s2mpjlib import *
class  UBH5(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : UBH5
#    *********
# 
#    The problem is to minimize the integral of the control magnitude needed
#    to bring a vehicle, from given position and velocity, to the origin with
#    zero velocity in a fixed amount of time.  The controls are the components
#    of the vehicle acceleration. The discretization uses the trapezoidal rule.
#    This version of the problem is a variant of UBH1, where the cumulative
#    value of the objective is maintained as an additional state variable.
# 
#    The problem is convex.
# 
#    Source: unscaled problem 5 
#    (ODE = 1, CLS = 2, GRD = 1, MET = T, SEED = 0.) in
#    J.T. Betts and W.P. Huffman,
#    "Sparse Nonlinear Programming Test Problems (Release 1.0)",
#    Boeing Computer services, Seattle, July 1993.
# 
#    SIF input: Ph.L. Toint, October 1993.
# 
#    classification = "LQR2-MN-V-V"
# 
#    Number of grid points
# 
#           Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER n=100, m=70    original value
# IE N                   100            $-PARAMETER n=1000, m=700
# IE N                   500            $-PARAMETER n=5000, m=3500
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'UBH5'

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
            v_['N'] = int(10);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
# IE N                   1000           $-PARAMETER n=10000, m=7000
# IE N                   2000           $-PARAMETER n=20000, m=14000
        v_['T0'] = 0.0
        v_['TF'] = 1000.0
        v_['RN'] = float(v_['N'])
        v_['TTIME'] = v_['TF']-v_['T0']
        v_['K'] = v_['TTIME']/v_['RN']
        v_['-K/2'] = -0.5*v_['K']
        v_['0'] = 0
        v_['1'] = 1
        v_['2'] = 2
        v_['3'] = 3
        v_['4'] = 4
        v_['5'] = 5
        v_['6'] = 6
        v_['7'] = 7
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        pb.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['7'])+1):
            for T in range(int(v_['0']),int(v_['N'])+1):
                [iv,ix_,_] = s2mpj_ii('Y'+str(I)+','+str(T),ix_)
                pb.xnames=arrset(pb.xnames,iv,'Y'+str(I)+','+str(T))
        for I in range(int(v_['1']),int(v_['3'])+1):
            for T in range(int(v_['0']),int(v_['N'])+1):
                [iv,ix_,_] = s2mpj_ii('U'+str(I)+','+str(T),ix_)
                pb.xnames=arrset(pb.xnames,iv,'U'+str(I)+','+str(T))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['Y'+str(int(v_['7']))+','+str(int(v_['N']))]
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['3'])+1):
            v_['I+3'] = 3+I
            for T in range(int(v_['1']),int(v_['N'])+1):
                v_['T-1'] = -1+T
                [ig,ig_,_] = s2mpj_ii('S'+str(I)+','+str(T),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'S'+str(I)+','+str(T))
                iv = ix_['Y'+str(I)+','+str(T)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
                iv = ix_['Y'+str(I)+','+str(int(v_['T-1']))]
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
                iv = ix_['Y'+str(int(v_['I+3']))+','+str(int(v_['T-1']))]
                pbm.A[ig,iv] = float(v_['-K/2'])+pbm.A[ig,iv]
                iv = ix_['Y'+str(int(v_['I+3']))+','+str(T)]
                pbm.A[ig,iv] = float(v_['-K/2'])+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['3'])+1):
            v_['I+3'] = 3+I
            for T in range(int(v_['1']),int(v_['N'])+1):
                v_['T-1'] = -1+T
                [ig,ig_,_] = s2mpj_ii('S'+str(int(v_['I+3']))+','+str(T),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'S'+str(int(v_['I+3']))+','+str(T))
                iv = ix_['Y'+str(int(v_['I+3']))+','+str(T)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
                iv = ix_['Y'+str(int(v_['I+3']))+','+str(int(v_['T-1']))]
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
                [ig,ig_,_] = s2mpj_ii('S'+str(int(v_['I+3']))+','+str(T),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'S'+str(int(v_['I+3']))+','+str(T))
                iv = ix_['U'+str(I)+','+str(int(v_['T-1']))]
                pbm.A[ig,iv] = float(v_['-K/2'])+pbm.A[ig,iv]
                [ig,ig_,_] = s2mpj_ii('S'+str(int(v_['I+3']))+','+str(T),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'S'+str(int(v_['I+3']))+','+str(T))
                iv = ix_['U'+str(I)+','+str(T)]
                pbm.A[ig,iv] = float(v_['-K/2'])+pbm.A[ig,iv]
        for T in range(int(v_['1']),int(v_['N'])+1):
            v_['T-1'] = -1+T
            [ig,ig_,_] = s2mpj_ii('S'+str(int(v_['7']))+','+str(T),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'S'+str(int(v_['7']))+','+str(T))
            iv = ix_['Y'+str(int(v_['7']))+','+str(T)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['Y'+str(int(v_['7']))+','+str(int(v_['T-1']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
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
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower = np.zeros((pb.n,1))
        for I in range(int(v_['1']),int(v_['3'])+1):
            for T in range(int(v_['0']),int(v_['N'])+1):
                pb.xlower[ix_['U'+str(I)+','+str(T)]] = -1.0
                pb.xupper[ix_['U'+str(I)+','+str(T)]] = 1.0
        pb.xlower[ix_['Y'+str(int(v_['1']))+','+str(int(v_['0']))]] = 1000.0
        pb.xupper[ix_['Y'+str(int(v_['1']))+','+str(int(v_['0']))]] = 1000.0
        pb.xlower[ix_['Y'+str(int(v_['2']))+','+str(int(v_['0']))]] = 1000.0
        pb.xupper[ix_['Y'+str(int(v_['2']))+','+str(int(v_['0']))]] = 1000.0
        pb.xlower[ix_['Y'+str(int(v_['3']))+','+str(int(v_['0']))]] = 1000.0
        pb.xupper[ix_['Y'+str(int(v_['3']))+','+str(int(v_['0']))]] = 1000.0
        pb.xlower[ix_['Y'+str(int(v_['4']))+','+str(int(v_['0']))]] = -10.0
        pb.xupper[ix_['Y'+str(int(v_['4']))+','+str(int(v_['0']))]] = -10.0
        pb.xlower[ix_['Y'+str(int(v_['5']))+','+str(int(v_['0']))]] = 10.0
        pb.xupper[ix_['Y'+str(int(v_['5']))+','+str(int(v_['0']))]] = 10.0
        pb.xlower[ix_['Y'+str(int(v_['6']))+','+str(int(v_['0']))]] = -10.0
        pb.xupper[ix_['Y'+str(int(v_['6']))+','+str(int(v_['0']))]] = -10.0
        pb.xlower[ix_['Y'+str(int(v_['7']))+','+str(int(v_['0']))]] = 0.0
        pb.xupper[ix_['Y'+str(int(v_['7']))+','+str(int(v_['0']))]] = 0.0
        for I in range(int(v_['1']),int(v_['6'])+1):
            pb.xlower[ix_['Y'+str(I)+','+str(int(v_['N']))]] = 0.0
            pb.xupper[ix_['Y'+str(I)+','+str(int(v_['N']))]] = 0.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'V')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for T in range(int(v_['0']),int(v_['N'])+1):
            for I in range(int(v_['1']),int(v_['3'])+1):
                ename = 'E'+str(I)+','+str(T)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
                ielftype = arrset(ielftype, ie, iet_["eSQ"])
                vname = 'U'+str(I)+','+str(T)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for T in range(int(v_['1']),int(v_['N'])+1):
            v_['T-1'] = -1+T
            for I in range(int(v_['1']),int(v_['3'])+1):
                ig = ig_['S'+str(int(v_['7']))+','+str(T)]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['E'+str(I)+','+str(int(v_['T-1']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-K/2']))
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(I)+','+str(T)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-K/2']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN(10)           1.14735202967
# LO SOLTN(100)          1.11631518169
# LO SOLTN(1000)         1.11598643493
# LO SOLTN(2000)         1.11587382445
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
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "LQR2-MN-V-V"
        pb.x0          = np.zeros((pb.n,1))
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

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

