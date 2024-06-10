from s2mpjlib import *
class  MADSSCHJ(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : MADSSCHJ
#    *********
# 
#    A nonlinear minmax problem with variable dimension.
#    The Jacobian of the constraints is dense.
# 
#    Source:
#    K. Madsen and H. Schjaer-Jacobsen,
#    "Linearly Constrained Minmax Optimization",
#    Mathematical Programming 14, pp. 208-223, 1978.
# 
#    SIF input: Ph. Toint, August 1993.
# 
#    classification = "LQR2-AN-V-V"
# 
#    N is the number of variables - 1, and must be even and at least 4.
#    The number of inequality constraints is 2*N - 2.
# 
#           Alternative values for the SIF file parameters:
# IE N                   4              $-PARAMETER  n=  5, m=  6
# IE N                   10             $-PARAMETER  n= 11, m= 18  original value
# IE N                   20             $-PARAMETER  n= 21, m= 38
# IE N                   30             $-PARAMETER  n= 31, m= 58
# IE N                   40             $-PARAMETER  n= 41, m= 78
# IE N                   50             $-PARAMETER  n= 51, m= 98
# IE N                   60             $-PARAMETER  n= 61, m=118
# IE N                   70             $-PARAMETER  n= 71, m=138
# IE N                   80             $-PARAMETER  n= 81, m=158
# IE N                   90             $-PARAMETER  n= 91, m=178
# IE N                   100            $-PARAMETER  n=101, m=198
# IE N                   200            $-PARAMETER  n=201, m=398
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'MADSSCHJ'

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
            v_['N'] = int(4);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
        v_['1'] = 1
        v_['2'] = 2
        v_['3'] = 3
        v_['4'] = 4
        v_['N-1'] = -1+v_['N']
        v_['2N'] = v_['N']+v_['N']
        v_['M'] = -2+v_['2N']
        v_['M-1'] = -1+v_['M']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
        [iv,ix_,_] = s2mpj_ii('Z',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Z')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['Z']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C1',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C1')
        iv = ix_['Z']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        for I in range(int(v_['2']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('C1',ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'C1')
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C2',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C2')
        iv = ix_['Z']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X1']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        for I in range(int(v_['3']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('C2',ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'C2')
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C3',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C3')
        iv = ix_['Z']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X1']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        for I in range(int(v_['3']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('C3',ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'C3')
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        for K in range(int(v_['4']),int(v_['M-1'])+1,int(v_['2'])):
            v_['K+1'] = 1+K
            v_['K+2'] = 2+K
            v_['J'] = int(np.fix(v_['K+2']/v_['2']))
            v_['J-1'] = -1+v_['J']
            v_['J+1'] = 1+v_['J']
            [ig,ig_,_] = s2mpj_ii('C'+str(K),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'C'+str(K))
            iv = ix_['Z']
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('C'+str(int(v_['K+1'])),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'C'+str(int(v_['K+1'])))
            iv = ix_['Z']
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            for I in range(int(v_['1']),int(v_['J-1'])+1):
                [ig,ig_,_] = s2mpj_ii('C'+str(K),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'C'+str(K))
                iv = ix_['X'+str(I)]
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
                [ig,ig_,_] = s2mpj_ii('C'+str(int(v_['K+1'])),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'C'+str(int(v_['K+1'])))
                iv = ix_['X'+str(I)]
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            for I in range(int(v_['J+1']),int(v_['N'])+1):
                [ig,ig_,_] = s2mpj_ii('C'+str(K),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'C'+str(K))
                iv = ix_['X'+str(I)]
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
                [ig,ig_,_] = s2mpj_ii('C'+str(int(v_['K+1'])),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'C'+str(int(v_['K+1'])))
                iv = ix_['X'+str(I)]
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C'+str(int(v_['M'])),ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C'+str(int(v_['M'])))
        iv = ix_['Z']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            [ig,ig_,_] = s2mpj_ii('C'+str(int(v_['M'])),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'C'+str(int(v_['M'])))
            iv = ix_['X'+str(I)]
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
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        for K in range(int(v_['1']),int(v_['M'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['C'+str(K)],float(-1.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower = np.zeros((pb.n,1))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(10.0))
        pb.x0[ix_['Z']] = float(0.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['1']),int(v_['N'])+1):
            ename = 'XSQ'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
            ielftype = arrset(ielftype, ie, iet_["eSQ"])
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,10.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        ig = ig_['C'+str(int(v_['1']))]
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['XSQ'+str(int(v_['1']))])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        ig = ig_['C'+str(int(v_['2']))]
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['XSQ'+str(int(v_['2']))])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        ig = ig_['C'+str(int(v_['3']))]
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['XSQ'+str(int(v_['2']))])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-2.0))
        for K in range(int(v_['4']),int(v_['M-1'])+1,int(v_['2'])):
            v_['K+1'] = 1+K
            v_['K+2'] = 2+K
            v_['J'] = int(np.fix(v_['K+2']/v_['2']))
            ig = ig_['C'+str(K)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['XSQ'+str(int(v_['J']))])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
            ig = ig_['C'+str(int(v_['K+1']))]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['XSQ'+str(int(v_['J']))])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(-2.0))
        ig = ig_['C'+str(int(v_['M']))]
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['XSQ'+str(int(v_['N']))])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN(4)             -2.6121094144
# LO SOLTN(10)            -12.814452425
# LO SOLTN(20)            -49.869888156
# LO SOLTN(30)            -111.93545559
# LO SOLTN(40)            -199.00371592
# LO SOLTN(50)            -311.07308068
# LO SOLTN(60)            -448.14300524
# LO SOLTN(70)            -610.21325256
# LO SOLTN(80)            -797.28370289
# LO SOLTN(90)            -1009.3542892
# LO SOLTN(100)           -1246.4249710
# LO SOLTN(200)           -4992.1339031
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle+pb.neq,pb.m)] = np.zeros((pb.nge,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "LQR2-AN-V-V"
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

