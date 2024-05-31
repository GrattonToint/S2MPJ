from s2xlib import *
class  KISSING(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem: KISSING NUMBER PROBLEM
#                                                                    
#    Source: This problem is associated to the family of Hard-Spheres 
#    problem. It belongs to the family of sphere packing problems, a 
#    class of challenging problems dating from the beginning of the 
#    17th century which is related to practical problems in Chemistry, 
#    Biology and Physics. It consists on maximizing the minimum pairwise 
#    distance between NP points on a sphere in \R^{MDIM}. 
#    This problem may be reduced to a nonconvex nonlinear optimization 
#    problem with a potentially large number of (nonoptimal) points 
#    satisfying optimality conditions. We have, thus, a class of problems 
#    indexed by the parameters MDIM and NP, that provides a suitable 
#    set of test problems for evaluating nonlinear programming codes.
#    After some algebric manipulations, we can formulate this problem as
#                             Minimize z
#                             subject to
#        
#       z \geq <x_i, x_j> for all different pair of indices i, j
#       
#                             ||x_i||^2 = 1    for all i = 1,...,NP
#      The goal is to find an objective value less than 0.5 (This means
#      that the NP points stored belong to the sphere and every distance
#      between two of them is greater than 1.0).
#      Obs: the starting point is aleatorally chosen although each 
#      variable belongs to [-1.,1.].
#      References:
#      [1] "Validation of an Augmented Lagrangian algorithm with a 
#           Gauss-Newton Hessian approximation using a set of 
#           Hard-Spheres problems", N. Krejic, J. M. Martinez, M. Mello 
#           and E. A. Pilotta, Tech. Report RP 29/98, IMECC-UNICAMP, 
#           Campinas, 1998.
#      [2] "Inexact-Restoration Algorithm for Constrained Optimization",
#           J. M. Martinez and E. A. Pilotta, Tech. Report, IMECC-UNICAMP, 
#           Campinas, 1998.
#      [3]  "Sphere Packings, Lattices and Groups", J. H. Conway and 
#            N. J. C. Sloane, Springer-Verlag, NY, 1988.
#      SIF input: September 29, 1998
# 		 Jose Mario Martinez
#                 Elvio Angel Pilotta
# 
#    classification = "LQR2-RN-V-V"
# 
# **********************************************************************
# 
#    Number of points: NP >= 12
# 
#           Alternative values for the SIF file parameters:
# IE NP                   12            $-PARAMETER
# IE NP                   13            $-PARAMETER
# IE NP                   14            $-PARAMETER
# IE NP                   15            $-PARAMETER
# IE NP                   22            $-PARAMETER
# IE NP                   23            $-PARAMETER
# IE NP                   24            $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'KISSING'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'KISSING'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['NP'] = int(25);  #  SIF file default value
        else:
            v_['NP'] = int(args[0])
# IE NP                   26            $-PARAMETER
# IE NP                   27            $-PARAMETER
# IE NP	                 37            $-PARAMETER
# IE NP                   38            $-PARAMETER
# IE NP                   39            $-PARAMETER
# IE NP                   40            $-PARAMETER
# IE NP                   41            $-PARAMETER
# IE NP                   42            $-PARAMETER
        if nargin<2:
            v_['MDIM'] = int(3);  #  SIF file default value
        else:
            v_['MDIM'] = int(args[1])
# IE MDIM                 4             $-PARAMETER
# IE MDIM                 5             $-PARAMETER
        v_['N-'] = -1+v_['NP']
        v_['1'] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['NP'])+1):
            for J in range(int(v_['1']),int(v_['MDIM'])+1):
                [iv,ix_,_] = s2x_ii('X'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'X'+str(I)+','+str(J))
        [iv,ix_,_] = s2x_ii('Z',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Z')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2x_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['Z']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['N-'])+1):
            v_['I+'] = 1+I
            for J in range(int(v_['I+']),int(v_['NP'])+1):
                [ig,ig_,_] = s2x_ii('IC'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<=')
                cnames = arrset(cnames,ig,'IC'+str(I)+','+str(J))
                iv = ix_['Z']
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['NP'])+1):
            [ig,ig_,_] = s2x_ii('EC'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'EC'+str(I))
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
        for I in range(int(v_['1']),int(v_['NP'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['EC'+str(I)],float(1.0))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('inf'))
        pb.xupper = np.full((pb.n,1),float('inf'))
        for I in range(int(v_['1']),int(v_['NP'])+1):
            for J in range(int(v_['1']),int(v_['MDIM'])+1):
                pb.xlower[ix_['X'+str(I)+','+str(J)]] = -float('Inf')
                pb.xupper[ix_['X'+str(I)+','+str(J)]] = +float('Inf')
        pb.xlower[ix_['Z']] = -float('Inf')
        pb.xupper[ix_['Z']] = +float('Inf')
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        pb.x0[ix_['X1,1']] = float(-0.10890604)
        pb.x0[ix_['X1,2']] = float(0.85395078)
        pb.x0[ix_['X1,3']] = float(-0.45461680)
        pb.x0[ix_['X2,1']] = float(0.49883922)
        pb.x0[ix_['X2,2']] = float(-0.18439316)
        pb.x0[ix_['X2,3']] = float(-0.04798594)
        pb.x0[ix_['X3,1']] = float(0.28262888)
        pb.x0[ix_['X3,2']] = float(-0.48054070)
        pb.x0[ix_['X3,3']] = float(0.46715332)
        pb.x0[ix_['X4,1']] = float(-0.00580106)
        pb.x0[ix_['X4,2']] = float(-0.49987584)
        pb.x0[ix_['X4,3']] = float(-0.44130302)
        pb.x0[ix_['X5,1']] = float(0.81712540)
        pb.x0[ix_['X5,2']] = float(-0.36874258)
        pb.x0[ix_['X5,3']] = float(-0.68321896)
        pb.x0[ix_['X6,1']] = float(0.29642426)
        pb.x0[ix_['X6,2']] = float(0.82315508)
        pb.x0[ix_['X6,3']] = float(0.35938150)
        pb.x0[ix_['X7,1']] = float(0.09215152)
        pb.x0[ix_['X7,2']] = float(-0.53564686)
        pb.x0[ix_['X7,3']] = float(0.00191436)
        pb.x0[ix_['X8,1']] = float(0.11700318)
        pb.x0[ix_['X8,2']] = float(0.96722760)
        pb.x0[ix_['X8,3']] = float(-0.14916438)
        pb.x0[ix_['X9,1']] = float(0.01791524)
        pb.x0[ix_['X9,2']] = float(0.17759446)
        pb.x0[ix_['X9,3']] = float(-0.61875872)
        pb.x0[ix_['X10,1']] = float(-0.63833630)
        pb.x0[ix_['X10,2']] = float(0.80830972)
        pb.x0[ix_['X10,3']] = float(0.45846734)
        pb.x0[ix_['X11,1']] = float(0.28446456)
        pb.x0[ix_['X11,2']] = float(0.45686938)
        pb.x0[ix_['X11,3']] = float(0.16368980)
        pb.x0[ix_['X12,1']] = float(0.76557382)
        pb.x0[ix_['X12,2']] = float(0.16700944)
        pb.x0[ix_['X12,3']] = float(-0.31647534)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'ePROD', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        [it,iet_,_] = s2x_ii( 'eQUA', iet_)
        elftv = loaset(elftv,it,0,'V')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['1']),int(v_['N-'])+1):
            v_['I+'] = 1+I
            for J in range(int(v_['I+']),int(v_['NP'])+1):
                for K in range(int(v_['1']),int(v_['MDIM'])+1):
                    ename = 'A'+str(I)+','+str(J)+','+str(K)
                    [ie,ie_,_] = s2x_ii(ename,ie_)
                    pbm.elftype = arrset(pbm.elftype,ie,'ePROD')
                    ielftype = arrset(ielftype, ie, iet_["ePROD"])
                    vname = 'X'+str(I)+','+str(K)
                    [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='X')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    vname = 'X'+str(J)+','+str(K)
                    [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['1']),int(v_['NP'])+1):
            for K in range(int(v_['1']),int(v_['MDIM'])+1):
                ename = 'B'+str(I)+','+str(K)
                [ie,ie_,_] = s2x_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eQUA')
                ielftype = arrset(ielftype, ie, iet_["eQUA"])
                vname = 'X'+str(I)+','+str(K)
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['N-'])+1):
            v_['I+'] = 1+I
            for J in range(int(v_['I+']),int(v_['NP'])+1):
                for K in range(int(v_['1']),int(v_['MDIM'])+1):
                    ig = ig_['IC'+str(I)+','+str(J)]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A'+str(I)+','+str(J)+','+str(K)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        for I in range(int(v_['1']),int(v_['NP'])+1):
            for K in range(int(v_['1']),int(v_['MDIM'])+1):
                ig = ig_['EC'+str(I)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['B'+str(I)+','+str(K)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.cupper[np.arange(pb.nle)] = np.zeros((pb.nle,1))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "LQR2-RN-V-V"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def ePROD(pbm,nargout,*args):

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
                H_[0,1] = 1.0
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eQUA(pbm,nargout,*args):

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
