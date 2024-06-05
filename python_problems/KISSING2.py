from s2mpjlib import *
class  KISSING2(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem: A second formulation of the KISSING NUMBER PROBLEM
#                                                                    
#    Source: This problem is associated to the family of Hard-Spheres 
#    problem. It belongs to the family of sphere packing problems, a 
#    class of challenging problems dating from the beginning of the 
#    17th century which is related to practical problems in Chemistry, 
#    Biology and Physics. Given a fixed unit sphere at the origin in R^n, 
#    the problem consists of arranging a further m unit spheres so that 
#    sum of the distances to these spheres is as small as possible.
#    This problem may be reduced to a nonconvex nonlinear optimization 
#    problem with a potentially large number of (nonoptimal) points 
#    satisfying optimality conditions. We have, thus, a class of problems 
#    indexed by the parameters m and n, that provides a suitable 
#    set of test problems for evaluating nonlinear programming codes.
#    After some algebric manipulations, we can formulate this problem as
#               m
#     Minimize sum <p_i,p_i> - m n
#              i=1
#     subject to
#        
#      <p_i - p_j, p_i - p_j> >= 4 for all different pair of indices i, j
#      and  
#      <p_i, p_i> >= 4 for all indices i
#   
#      as well as n(n-1)/2 normalisation constraints fixing components.
#      The goal is to find an objective value equal to 0.
#      [1]  "Sphere Packings, Lattices and Groups", J. H. Conway and 
#            N. J. C. Sloane, Springer-Verlag, NY, 1988.
#    SIF input: Nick Gould, September 2000
# 
#    classification = "QQR2-RN-V-V"
# 
# **********************************************************************
# 
#    Number of points: m
# 
#           Alternative values for the SIF file parameters:
# IE m                   24             $-PARAMETER  number of points
# IE m                   25             $-PARAMETER  number of points
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'KISSING2'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'KISSING2'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['m'] = int(25);  #  SIF file default value
        else:
            v_['m'] = int(args[0])
# IE m                   100            $-PARAMETER  number of points
# IE n                    4             $-PARAMETER  dimension of sphere
        if nargin<2:
            v_['n'] = int(4);  #  SIF file default value
        else:
            v_['n'] = int(args[1])
# IE n                    8             $-PARAMETER  dimension of sphere
        v_['1'] = 1
        v_['2'] = 2
        v_['n-1'] = v_['n']-v_['1']
        v_['rm'] = float(v_['m'])
        v_['rn'] = float(v_['n'])
        v_['RM+N'] = v_['rm']+v_['rn']
        v_['mn'] = v_['rm']*v_['rn']
        v_['PI/4'] = np.arctan(1.0)
        v_['PI'] = 4.0*v_['PI/4']
        v_['PI/m'] = v_['PI']/v_['rm']
        v_['2PI/m'] = 2.0*v_['PI/m']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['m'])+1):
            for J in range(int(v_['1']),int(v_['n'])+1):
                [iv,ix_,_] = s2mpj_ii('P'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'P'+str(I)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['1']),int(v_['m'])+1):
            for J in range(int(v_['1']),int(v_['m'])+1):
                [ig,ig_,_] = s2mpj_ii('C'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'C'+str(I)+','+str(J))
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
        pbm.gconst = arrset(pbm.gconst,ig_['OBJ'],float(v_['mn']))
        for I in range(int(v_['1']),int(v_['m'])+1):
            for J in range(int(v_['1']),int(v_['m'])+1):
                pbm.gconst = arrset(pbm.gconst,ig_['C'+str(I)+','+str(J)],float(4.0))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),float('inf'))
        for I in range(int(v_['1']),int(v_['m'])+1):
            for J in range(int(v_['1']),int(v_['n'])+1):
                pb.xlower[ix_['P'+str(I)+','+str(J)]] = -float('Inf')
                pb.xupper[ix_['P'+str(I)+','+str(J)]] = +float('Inf')
        for I in range(int(v_['2']),int(v_['n'])+1):
            for J in range(int(I),int(v_['n'])+1):
                pb.xlower[ix_['P'+str(I)+','+str(J)]] = 0.0
                pb.xupper[ix_['P'+str(I)+','+str(J)]] = 0.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        for I in range(int(v_['1']),int(v_['m'])+1):
            v_['RI'] = float(I)
            v_['2PIi/m'] = v_['2PI/m']*v_['RI']
            v_['cos'] = np.cos(v_['2PIi/m'])
            v_['sin'] = np.sin(v_['2PIi/m'])
            v_['cos'] = v_['cos']+v_['cos']
            v_['sin'] = v_['sin']+v_['sin']
            if('P'+str(I)+','+str(int(v_['1'])) in ix_):
                pb.x0[ix_['P'+str(I)+','+str(int(v_['1']))]] = float(v_['cos'])
            else:
                pb.y0  = (
                      arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['P'+str(I)+','+str(int(v_['1']))]),float(v_['cos'])))
            for J in range(int(v_['2']),int(v_['n-1'])+1):
                if('P'+str(I)+','+str(J) in ix_):
                    pb.x0[ix_['P'+str(I)+','+str(J)]] = float(v_['sin'])
                else:
                    pb.y0  = (
                          arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['P'+str(I)+','+str(J)]),float(v_['sin'])))
        if('P'+str(int(v_['m']))+','+str(int(v_['n'])) in ix_):
            pb.x0[ix_['P'+str(int(v_['m']))+','+str(int(v_['n']))]] = float(v_['cos'])
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['P'+str(int(v_['m']))+','+str(int(v_['n']))]),float(v_['cos'])))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'ePROD1', iet_)
        elftv = loaset(elftv,it,0,'P')
        [it,iet_,_] = s2mpj_ii( 'ePROD2', iet_)
        elftv = loaset(elftv,it,0,'Q')
        elftv = loaset(elftv,it,1,'R')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['1']),int(v_['m'])+1):
            v_['I-'] = -1+I
            v_['I+'] = 1+I
            for J in range(int(v_['1']),int(v_['I-'])+1):
                for K in range(int(v_['1']),int(v_['n'])+1):
                    ename = 'E'+str(I)+','+str(J)+','+str(K)
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    pbm.elftype = arrset(pbm.elftype,ie,'ePROD2')
                    ielftype = arrset(ielftype, ie, iet_["ePROD2"])
                    vname = 'P'+str(I)+','+str(K)
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='Q')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    vname = 'P'+str(J)+','+str(K)
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='R')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            for K in range(int(v_['1']),int(v_['n'])+1):
                ename = 'E'+str(I)+','+str(I)+','+str(K)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'ePROD1')
                ielftype = arrset(ielftype, ie, iet_["ePROD1"])
                vname = 'P'+str(I)+','+str(K)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='P')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            for J in range(int(v_['I+']),int(v_['m'])+1):
                for K in range(int(v_['1']),int(v_['n'])+1):
                    ename = 'E'+str(I)+','+str(J)+','+str(K)
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    pbm.elftype = arrset(pbm.elftype,ie,'ePROD2')
                    ielftype = arrset(ielftype, ie, iet_["ePROD2"])
                    vname = 'P'+str(I)+','+str(K)
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='Q')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    vname = 'P'+str(J)+','+str(K)
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='R')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['m'])+1):
            for K in range(int(v_['1']),int(v_['n'])+1):
                ig = ig_['OBJ']
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(I)+','+str(I)+','+str(K)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            for J in range(int(v_['1']),int(v_['m'])+1):
                for K in range(int(v_['1']),int(v_['n'])+1):
                    ig = ig_['C'+str(I)+','+str(J)]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['E'+str(I)+','+str(J)+','+str(K)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# XL SOLUTION             0.00000D+00   $ n=4, m = 24
# XL SOLUTION             6.48030D+00   $ n=4, m = 25 one of many local solutions
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle+pb.neq,pb.m)] = np.zeros((pb.nge,1))
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "QQR2-RN-V-V"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def ePROD1(pbm,nargout,*args):

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

    @staticmethod
    def ePROD2(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((1,2))
        IV_ = np.zeros(1)
        U_[0,0] = U_[0,0]+1
        U_[0,1] = U_[0,1]-1
        IV_[0] = U_[0:1,:].dot(EV_)
        f_   = IV_[0]*IV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = IV_[0]+IV_[0]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

