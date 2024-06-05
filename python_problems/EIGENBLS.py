from s2mpjlib import *
class  EIGENBLS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : EIGENBLS
#    --------
# 
#    Solving symmetric eigenvalue problems as systems of
#    nonlinear equations.
# 
#    The problem is, given a symmetric matrix A, to find an orthogonal
#    matrix Q and diagonal matrix D such that A = Q(T) D Q.
# 
#    Example B: a tridiagonal matrix with diagonals 2 and off diagonals -1
# 
#    Source:  An idea by Nick Gould
# 
#             Least-squares version
# 
#    SIF input: Nick Gould, Nov 1992.
# 
#    classification = "SUR2-AN-V-0"
# 
#    The dimension of the matrix.
# 
#           Alternative values for the SIF file parameters:
# IE N                   2              $-PARAMETER
# IE N                   10             $-PARAMETER     original value
# IE N                   50             $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'EIGENBLS'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'EIGENBLS'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(2);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
        v_['1'] = 1
        v_['2'] = 2
        v_['N-1'] = -1+v_['N']
        v_['A'+str(int(v_['1']))+','+str(int(v_['1']))] = 2.0
        for J in range(int(v_['2']),int(v_['N'])+1):
            v_['J-1'] = -1+J
            v_['J-2'] = -2+J
            for I in range(int(v_['1']),int(v_['J-2'])+1):
                v_['A'+str(I)+','+str(J)] = 0.0
            v_['A'+str(int(v_['J-1']))+','+str(J)] = -1.0
            v_['A'+str(J)+','+str(J)] = 2.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for J in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('D'+str(J),ix_)
            pb.xnames=arrset(pb.xnames,iv,'D'+str(J))
            for I in range(int(v_['1']),int(v_['N'])+1):
                [iv,ix_,_] = s2mpj_ii('Q'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'Q'+str(I)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for J in range(int(v_['1']),int(v_['N'])+1):
            for I in range(int(v_['1']),int(J)+1):
                [ig,ig_,_] = s2mpj_ii('E'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
                [ig,ig_,_] = s2mpj_ii('O'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        for J in range(int(v_['1']),int(v_['N'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['O'+str(J)+','+str(J)],float(1.0))
            for I in range(int(v_['1']),int(J)+1):
                pbm.gconst  = (
                      arrset(pbm.gconst,ig_['E'+str(I)+','+str(J)],float(v_['A'+str(I)+','+str(J)])))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(0.0))
        for J in range(int(v_['1']),int(v_['N'])+1):
            pb.x0[ix_['D'+str(J)]] = float(1.0)
            pb.x0[ix_['Q'+str(J)+','+str(J)]] = float(1.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'en2PROD', iet_)
        elftv = loaset(elftv,it,0,'Q1')
        elftv = loaset(elftv,it,1,'Q2')
        [it,iet_,_] = s2mpj_ii( 'en3PROD', iet_)
        elftv = loaset(elftv,it,0,'Q1')
        elftv = loaset(elftv,it,1,'Q2')
        elftv = loaset(elftv,it,2,'D')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for J in range(int(v_['1']),int(v_['N'])+1):
            for I in range(int(v_['1']),int(J)+1):
                for K in range(int(v_['1']),int(v_['N'])+1):
                    ename = 'E'+str(I)+','+str(J)+','+str(K)
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    pbm.elftype = arrset(pbm.elftype,ie,'en3PROD')
                    ielftype = arrset(ielftype, ie, iet_["en3PROD"])
                    vname = 'Q'+str(K)+','+str(I)
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='Q1')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    vname = 'Q'+str(K)+','+str(J)
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='Q2')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    vname = 'D'+str(K)
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='D')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    ename = 'O'+str(I)+','+str(J)+','+str(K)
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    pbm.elftype = arrset(pbm.elftype,ie,'en2PROD')
                    ielftype = arrset(ielftype, ie, iet_["en2PROD"])
                    vname = 'Q'+str(K)+','+str(I)
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='Q1')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    vname = 'Q'+str(K)+','+str(J)
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='Q2')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
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
        for ig in range(0,ngrp):
            pbm.grftype = arrset(pbm.grftype,ig,'gL2')
        for J in range(int(v_['1']),int(v_['N'])+1):
            for I in range(int(v_['1']),int(J)+1):
                for K in range(int(v_['1']),int(v_['N'])+1):
                    ig = ig_['E'+str(I)+','+str(J)]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['E'+str(I)+','+str(J)+','+str(K)]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
                    ig = ig_['O'+str(I)+','+str(J)]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['O'+str(I)+','+str(J)+','+str(K)]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "SUR2-AN-V-0"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def en2PROD(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]
            g_[1] = EV_[0]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = 1.0e+0
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def en3PROD(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[1]*EV_[2]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]*EV_[2]
            g_[1] = EV_[0]*EV_[2]
            g_[2] = EV_[0]*EV_[1]
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = EV_[2]
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1]
                H_[2,0] = H_[0,2]
                H_[1,2] = EV_[0]
                H_[2,1] = H_[1,2]
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
                H_ = 2.0e+0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

