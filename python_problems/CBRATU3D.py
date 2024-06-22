from s2mpjlib import *
class  CBRATU3D(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CBRATU3D
#    *********
# 
#    The complex 3D Bratu problem on the unit cube, using finite
#    differences.
# 
#    Source: Problem 3 in
#    J.J. More',
#    "A collection of nonlinear model problems"
#    Proceedings of the AMS-SIAM Summer seminar on the Computational
#    Solution of Nonlinear Systems of Equations, Colorado, 1988.
#    Argonne National Laboratory MCS-P60-0289, 1989.
# 
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "NOR2-MN-V-V"
# 
#    P is the number of points in one side of the unit cube
#    There are 2*P**3 variables
# 
#           Alternative values for the SIF file parameters:
# IE P                   3              $-PARAMETER n = 54   original value
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'CBRATU3D'

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
            v_['P'] = int(3);  #  SIF file default value
        else:
            v_['P'] = int(args[0])
# IE P                   4              $-PARAMETER n = 128
# IE P                   7              $-PARAMETER n = 686
# IE P                   10             $-PARAMETER n = 2000
# IE P                   12             $-PARAMETER n = 3456
        if nargin<2:
            v_['LAMBDA'] = float(6.80812);  #  SIF file default value
        else:
            v_['LAMBDA'] = float(args[1])
        v_['1'] = 1
        v_['2'] = 2
        v_['1.0'] = 1.0
        v_['P-1'] = -1+v_['P']
        v_['RP-1'] = float(v_['P-1'])
        v_['H'] = v_['1.0']/v_['RP-1']
        v_['H2'] = v_['H']*v_['H']
        v_['C'] = v_['H2']*v_['LAMBDA']
        v_['-C'] = -1.0*v_['C']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        pb.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for J in range(int(v_['1']),int(v_['P'])+1):
            for I in range(int(v_['1']),int(v_['P'])+1):
                for K in range(int(v_['1']),int(v_['P'])+1):
                    [iv,ix_,_] = s2mpj_ii('U'+str(I)+','+str(J)+','+str(K),ix_)
                    pb.xnames=arrset(pb.xnames,iv,'U'+str(I)+','+str(J)+','+str(K))
                    [iv,ix_,_] = s2mpj_ii('X'+str(I)+','+str(J)+','+str(K),ix_)
                    pb.xnames=arrset(pb.xnames,iv,'X'+str(I)+','+str(J)+','+str(K))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['2']),int(v_['P-1'])+1):
            v_['R'] = 1+I
            v_['S'] = -1+I
            for J in range(int(v_['2']),int(v_['P-1'])+1):
                v_['V'] = 1+J
                v_['W'] = -1+J
                for K in range(int(v_['2']),int(v_['P-1'])+1):
                    v_['Y'] = 1+K
                    v_['Z'] = -1+K
                    [ig,ig_,_] = s2mpj_ii('G'+str(I)+','+str(J)+','+str(K),ig_)
                    gtype = arrset(gtype,ig,'==')
                    cnames = arrset(cnames,ig,'G'+str(I)+','+str(J)+','+str(K))
                    iv = ix_['U'+str(I)+','+str(J)+','+str(K)]
                    pbm.A[ig,iv] = float(6.0)+pbm.A[ig,iv]
                    iv = ix_['U'+str(int(v_['R']))+','+str(J)+','+str(K)]
                    pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
                    iv = ix_['U'+str(int(v_['S']))+','+str(J)+','+str(K)]
                    pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
                    iv = ix_['U'+str(I)+','+str(int(v_['V']))+','+str(K)]
                    pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
                    iv = ix_['U'+str(I)+','+str(int(v_['W']))+','+str(K)]
                    pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
                    iv = ix_['U'+str(I)+','+str(J)+','+str(int(v_['Y']))]
                    pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
                    iv = ix_['U'+str(I)+','+str(J)+','+str(int(v_['Z']))]
                    pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
                    [ig,ig_,_] = s2mpj_ii('F'+str(I)+','+str(J)+','+str(K),ig_)
                    gtype = arrset(gtype,ig,'==')
                    cnames = arrset(cnames,ig,'F'+str(I)+','+str(J)+','+str(K))
                    iv = ix_['X'+str(I)+','+str(J)+','+str(K)]
                    pbm.A[ig,iv] = float(6.0)+pbm.A[ig,iv]
                    iv = ix_['X'+str(int(v_['R']))+','+str(J)+','+str(K)]
                    pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
                    iv = ix_['X'+str(int(v_['S']))+','+str(J)+','+str(K)]
                    pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
                    iv = ix_['X'+str(I)+','+str(int(v_['V']))+','+str(K)]
                    pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
                    iv = ix_['X'+str(I)+','+str(int(v_['W']))+','+str(K)]
                    pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
                    iv = ix_['X'+str(I)+','+str(J)+','+str(int(v_['Y']))]
                    pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
                    iv = ix_['X'+str(I)+','+str(J)+','+str(int(v_['Z']))]
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
        for J in range(int(v_['1']),int(v_['P'])+1):
            for K in range(int(v_['1']),int(v_['P'])+1):
                pb.xlower[ix_['U'+str(int(v_['1']))+','+str(J)+','+str(K)]] = 0.0
                pb.xupper[ix_['U'+str(int(v_['1']))+','+str(J)+','+str(K)]] = 0.0
                pb.xlower[ix_['U'+str(int(v_['P']))+','+str(J)+','+str(K)]] = 0.0
                pb.xupper[ix_['U'+str(int(v_['P']))+','+str(J)+','+str(K)]] = 0.0
                pb.xlower[ix_['X'+str(int(v_['1']))+','+str(J)+','+str(K)]] = 0.0
                pb.xupper[ix_['X'+str(int(v_['1']))+','+str(J)+','+str(K)]] = 0.0
                pb.xlower[ix_['X'+str(int(v_['P']))+','+str(J)+','+str(K)]] = 0.0
                pb.xupper[ix_['X'+str(int(v_['P']))+','+str(J)+','+str(K)]] = 0.0
        for I in range(int(v_['2']),int(v_['P-1'])+1):
            for K in range(int(v_['1']),int(v_['P'])+1):
                pb.xlower[ix_['U'+str(I)+','+str(int(v_['P']))+','+str(K)]] = 0.0
                pb.xupper[ix_['U'+str(I)+','+str(int(v_['P']))+','+str(K)]] = 0.0
                pb.xlower[ix_['U'+str(I)+','+str(int(v_['1']))+','+str(K)]] = 0.0
                pb.xupper[ix_['U'+str(I)+','+str(int(v_['1']))+','+str(K)]] = 0.0
                pb.xlower[ix_['X'+str(I)+','+str(int(v_['P']))+','+str(K)]] = 0.0
                pb.xupper[ix_['X'+str(I)+','+str(int(v_['P']))+','+str(K)]] = 0.0
                pb.xlower[ix_['X'+str(I)+','+str(int(v_['1']))+','+str(K)]] = 0.0
                pb.xupper[ix_['X'+str(I)+','+str(int(v_['1']))+','+str(K)]] = 0.0
        for I in range(int(v_['2']),int(v_['P-1'])+1):
            for J in range(int(v_['2']),int(v_['P-1'])+1):
                pb.xlower[ix_['U'+str(I)+','+str(J)+','+str(int(v_['1']))]] = 0.0
                pb.xupper[ix_['U'+str(I)+','+str(J)+','+str(int(v_['1']))]] = 0.0
                pb.xlower[ix_['U'+str(I)+','+str(J)+','+str(int(v_['P']))]] = 0.0
                pb.xupper[ix_['U'+str(I)+','+str(J)+','+str(int(v_['P']))]] = 0.0
                pb.xlower[ix_['X'+str(I)+','+str(J)+','+str(int(v_['1']))]] = 0.0
                pb.xupper[ix_['X'+str(I)+','+str(J)+','+str(int(v_['1']))]] = 0.0
                pb.xlower[ix_['X'+str(I)+','+str(J)+','+str(int(v_['P']))]] = 0.0
                pb.xupper[ix_['X'+str(I)+','+str(J)+','+str(int(v_['P']))]] = 0.0
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(0.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eRPART', iet_)
        elftv = loaset(elftv,it,0,'U')
        elftv = loaset(elftv,it,1,'V')
        [it,iet_,_] = s2mpj_ii( 'eCPART', iet_)
        elftv = loaset(elftv,it,0,'U')
        elftv = loaset(elftv,it,1,'V')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['2']),int(v_['P-1'])+1):
            for J in range(int(v_['2']),int(v_['P-1'])+1):
                for K in range(int(v_['2']),int(v_['P-1'])+1):
                    ename = 'A'+str(I)+','+str(J)+','+str(K)
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    pbm.elftype = arrset(pbm.elftype,ie,'eRPART')
                    ielftype = arrset(ielftype, ie, iet_["eRPART"])
                    vname = 'U'+str(I)+','+str(J)+','+str(K)
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='U')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    vname = 'X'+str(I)+','+str(J)+','+str(K)
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='V')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    ename = 'B'+str(I)+','+str(J)+','+str(K)
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    pbm.elftype = arrset(pbm.elftype,ie,'eCPART')
                    ielftype = arrset(ielftype, ie, iet_["eCPART"])
                    vname = 'U'+str(I)+','+str(J)+','+str(K)
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='U')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    vname = 'X'+str(I)+','+str(J)+','+str(K)
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='V')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['2']),int(v_['P-1'])+1):
            for J in range(int(v_['2']),int(v_['P-1'])+1):
                for K in range(int(v_['2']),int(v_['P-1'])+1):
                    ig = ig_['G'+str(I)+','+str(J)+','+str(K)]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A'+str(I)+','+str(J)+','+str(K)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-C']))
                    ig = ig_['F'+str(I)+','+str(J)+','+str(K)]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B'+str(I)+','+str(J)+','+str(K)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-C']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN               0.0
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
        pb.pbclass = "NOR2-MN-V-V"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eRPART(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        EXPU = np.exp(EV_[0])
        EXPUC = EXPU*np.cos(EV_[1])
        EXPUS = EXPU*np.sin(EV_[1])
        f_   = EXPUC
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EXPUC
            g_[1] = -EXPUS
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = EXPUC
                H_[0,1] = -EXPUS
                H_[1,0] = H_[0,1]
                H_[1,1] = -EXPUC
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eCPART(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        EXPU = np.exp(EV_[0])
        EXPUC = EXPU*np.cos(EV_[1])
        EXPUS = EXPU*np.sin(EV_[1])
        f_   = EXPUS
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EXPUS
            g_[1] = EXPUC
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = EXPUS
                H_[0,1] = EXPUC
                H_[1,0] = H_[0,1]
                H_[1,1] = -EXPUS
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

