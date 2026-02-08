from s2mpjlib import *
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
#    classification = "C-CLQR2-RN-V-V"
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
# IE NP                   25            $-PARAMETER
# IE NP                   26            $-PARAMETER
# IE NP                   27            $-PARAMETER
# IE NP	                 37            $-PARAMETER
# IE NP                   38            $-PARAMETER
# IE NP                   39            $-PARAMETER
# IE NP                   40            $-PARAMETER
# IE NP                   41            $-PARAMETER
# IE NP                   42            $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'KISSING'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['NP'] = int(25);  #  SIF file default value
        else:
            v_['NP'] = int(args[0])
# IE MDIM                 3             $-PARAMETER
        if nargin<2:
            v_['MDIM'] = int(3);  #  SIF file default value
        else:
            v_['MDIM'] = int(args[1])
# IE MDIM                 4             $-PARAMETER
# IE MDIM                 5             $-PARAMETER
        v_['N-'] = -1+v_['NP']
        v_['1'] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['NP'])+1):
            for J in range(int(v_['1']),int(v_['MDIM'])+1):
                [iv,ix_,_] = s2mpj_ii('X'+str(I)+','+str(J),ix_)
                self.xnames=arrset(self.xnames,iv,'X'+str(I)+','+str(J))
        [iv,ix_,_] = s2mpj_ii('Z',ix_)
        self.xnames=arrset(self.xnames,iv,'Z')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['Z']])
        valA = np.append(valA,float(1.0))
        for I in range(int(v_['1']),int(v_['N-'])+1):
            v_['I+'] = 1+I
            for J in range(int(v_['I+']),int(v_['NP'])+1):
                [ig,ig_,_] = s2mpj_ii('IC'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<=')
                cnames = arrset(cnames,ig,'IC'+str(I)+','+str(J))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['Z']])
                valA = np.append(valA,float(-1.0))
        for I in range(int(v_['1']),int(v_['NP'])+1):
            [ig,ig_,_] = s2mpj_ii('EC'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'EC'+str(I))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        legrps = np.where(gtype=='<=')[0]
        eqgrps = np.where(gtype=='==')[0]
        gegrps = np.where(gtype=='>=')[0]
        self.nle = len(legrps)
        self.neq = len(eqgrps)
        self.nge = len(gegrps)
        self.m   = self.nle+self.neq+self.nge
        self.congrps = np.concatenate((legrps,eqgrps,gegrps))
        self.cnames = cnames[self.congrps]
        self.nob = ngrp-self.m
        self.objgrps = np.where(gtype=='<>')[0]
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        for I in range(int(v_['1']),int(v_['NP'])+1):
            self.gconst = arrset(self.gconst,ig_['EC'+str(I)],float(1.0))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        for I in range(int(v_['1']),int(v_['NP'])+1):
            for J in range(int(v_['1']),int(v_['MDIM'])+1):
                self.xlower[ix_['X'+str(I)+','+str(J)]] = -float('Inf')
                self.xupper[ix_['X'+str(I)+','+str(J)]] = +float('Inf')
        self.xlower[ix_['Z']] = -float('Inf')
        self.xupper[ix_['Z']] = +float('Inf')
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        self.x0[ix_['X1,1']] = float(-0.10890604)
        self.x0[ix_['X1,2']] = float(0.85395078)
        self.x0[ix_['X1,3']] = float(-0.45461680)
        self.x0[ix_['X2,1']] = float(0.49883922)
        self.x0[ix_['X2,2']] = float(-0.18439316)
        self.x0[ix_['X2,3']] = float(-0.04798594)
        self.x0[ix_['X3,1']] = float(0.28262888)
        self.x0[ix_['X3,2']] = float(-0.48054070)
        self.x0[ix_['X3,3']] = float(0.46715332)
        self.x0[ix_['X4,1']] = float(-0.00580106)
        self.x0[ix_['X4,2']] = float(-0.49987584)
        self.x0[ix_['X4,3']] = float(-0.44130302)
        self.x0[ix_['X5,1']] = float(0.81712540)
        self.x0[ix_['X5,2']] = float(-0.36874258)
        self.x0[ix_['X5,3']] = float(-0.68321896)
        self.x0[ix_['X6,1']] = float(0.29642426)
        self.x0[ix_['X6,2']] = float(0.82315508)
        self.x0[ix_['X6,3']] = float(0.35938150)
        self.x0[ix_['X7,1']] = float(0.09215152)
        self.x0[ix_['X7,2']] = float(-0.53564686)
        self.x0[ix_['X7,3']] = float(0.00191436)
        self.x0[ix_['X8,1']] = float(0.11700318)
        self.x0[ix_['X8,2']] = float(0.96722760)
        self.x0[ix_['X8,3']] = float(-0.14916438)
        self.x0[ix_['X9,1']] = float(0.01791524)
        self.x0[ix_['X9,2']] = float(0.17759446)
        self.x0[ix_['X9,3']] = float(-0.61875872)
        self.x0[ix_['X10,1']] = float(-0.63833630)
        self.x0[ix_['X10,2']] = float(0.80830972)
        self.x0[ix_['X10,3']] = float(0.45846734)
        self.x0[ix_['X11,1']] = float(0.28446456)
        self.x0[ix_['X11,2']] = float(0.45686938)
        self.x0[ix_['X11,3']] = float(0.16368980)
        self.x0[ix_['X12,1']] = float(0.76557382)
        self.x0[ix_['X12,2']] = float(0.16700944)
        self.x0[ix_['X12,3']] = float(-0.31647534)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'ePROD', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        [it,iet_,_] = s2mpj_ii( 'eQUA', iet_)
        elftv = loaset(elftv,it,0,'V')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['1']),int(v_['N-'])+1):
            v_['I+'] = 1+I
            for J in range(int(v_['I+']),int(v_['NP'])+1):
                for K in range(int(v_['1']),int(v_['MDIM'])+1):
                    ename = 'A'+str(I)+','+str(J)+','+str(K)
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    self.elftype = arrset(self.elftype,ie,'ePROD')
                    ielftype = arrset(ielftype,ie,iet_["ePROD"])
                    vname = 'X'+str(I)+','+str(K)
                    [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                    posev = np.where(elftv[ielftype[ie]]=='X')[0]
                    self.elvar = loaset(self.elvar,ie,posev[0],iv)
                    vname = 'X'+str(J)+','+str(K)
                    [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                    posev = np.where(elftv[ielftype[ie]]=='Y')[0]
                    self.elvar = loaset(self.elvar,ie,posev[0],iv)
        for I in range(int(v_['1']),int(v_['NP'])+1):
            for K in range(int(v_['1']),int(v_['MDIM'])+1):
                ename = 'B'+str(I)+','+str(K)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eQUA')
                ielftype = arrset(ielftype,ie,iet_["eQUA"])
                vname = 'X'+str(I)+','+str(K)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                posev = np.where(elftv[ielftype[ie]]=='V')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['N-'])+1):
            v_['I+'] = 1+I
            for J in range(int(v_['I+']),int(v_['NP'])+1):
                for K in range(int(v_['1']),int(v_['MDIM'])+1):
                    ig = ig_['IC'+str(I)+','+str(J)]
                    posel = len(self.grelt[ig])
                    self.grelt  = (
                          loaset(self.grelt,ig,posel,ie_['A'+str(I)+','+str(J)+','+str(K)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    self.grelw = loaset(self.grelw,ig,posel,1.)
        for I in range(int(v_['1']),int(v_['NP'])+1):
            for K in range(int(v_['1']),int(v_['MDIM'])+1):
                ig = ig_['EC'+str(I)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['B'+str(I)+','+str(K)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# XL SOLUTION             4.47214D-01
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.cupper[np.arange(self.nle)] = np.zeros((self.nle,1))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CLQR2-RN-V-V"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def ePROD(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0,0]*EV_[1,0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1,0]
            g_[1] = EV_[0,0]
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
    def eQUA(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0,0]*EV_[0,0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[0,0]+EV_[0,0]
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

