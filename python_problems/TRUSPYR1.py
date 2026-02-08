from s2mpjlib import *
class  TRUSPYR1(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    This is a structural optimization problem.
#    The problem is to minimize the weight of a given
#    8-bar truss structure formed as a pyramid for a given external load.
#    There are upper bounds on the strain energy and lower bounds
#    on the cross-sectional areas of the bars.
# 
#    Source:
#    K. Svanberg, 
#    "Local and global optima",
#    Proceedings of the NATO/DFG ASI on Optimization of large structural
#    systems, 
#    G. I. N. Rozvany, ed., Kluwer, 1993, pp. 579-588.
# 
#    SIF input: A. Forsgren, Royal Institute of Technology, December 1993.
# 
#    classification = "C-CLQR2-MN-11-4"
# 
#    Number of bars
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'TRUSPYR1'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['NBAR'] = 8
        v_['NDIM'] = 3
        v_['1'] = 1
        v_['2'] = 2
        v_['NBAR/2'] = int(np.fix(v_['NBAR']/v_['2']))
        v_['8.0'] = 8.0
        v_['SQRT17'] = np.sqrt(17.0)
        v_['SQRT18'] = np.sqrt(18.0)
        v_['SQRT105'] = np.sqrt(105.0)
        v_['P1'] = 40.0
        v_['P2'] = 20.0
        v_['P3'] = 200.0
        v_['Q1'] = 2.0/v_['SQRT105']
        v_['Q2'] = 1.0/v_['SQRT105']
        v_['Q3'] = 10.0/v_['SQRT105']
        v_['ALPHA'] = 0.3291726437
        for J in range(int(v_['1']),int(v_['NBAR/2'])+1):
            v_['L'+str(J)] = v_['SQRT17']/v_['8.0']
            v_['J+4'] = J+v_['NBAR/2']
            v_['L'+str(int(v_['J+4']))] = v_['SQRT18']/v_['8.0']
        v_['E'] = 21.0
        v_['R1,1'] = 0.250
        v_['R2,1'] = 0.250
        v_['R3,1'] = 0.375
        v_['R1,2'] = 0.250
        v_['R2,2'] = -0.250
        v_['R3,2'] = 0.375
        v_['R1,3'] = -0.250
        v_['R2,3'] = -0.250
        v_['R3,3'] = 0.375
        v_['R1,4'] = -0.250
        v_['R2,4'] = 0.250
        v_['R3,4'] = 0.375
        v_['R1,5'] = 0.375
        v_['R2,5'] = 0.000
        v_['R3,5'] = 0.375
        v_['R1,6'] = 0.000
        v_['R2,6'] = -0.375
        v_['R3,6'] = 0.375
        v_['R1,7'] = -0.375
        v_['R2,7'] = 0.000
        v_['R3,7'] = 0.375
        v_['R1,8'] = 0.000
        v_['R2,8'] = 0.375
        v_['R3,8'] = 0.375
        for J in range(int(v_['1']),int(v_['NBAR'])+1):
            v_['L2'+str(J)] = v_['L'+str(J)]*v_['L'+str(J)]
            v_['L3'+str(J)] = v_['L2'+str(J)]*v_['L'+str(J)]
            v_['GAMMA'+str(J)] = v_['E']/v_['L3'+str(J)]
            v_['DL2'+str(J)] = v_['L2'+str(J)]/v_['E']
            v_['W'+str(J)] = 0.78*v_['L'+str(J)]
            v_['STRUP'+str(J)] = 10.0*v_['DL2'+str(J)]
            for I in range(int(v_['1']),int(v_['NDIM'])+1):
                v_['RG'+str(I)+','+str(J)] = v_['GAMMA'+str(J)]*v_['R'+str(I)+','+str(J)]
                for K in range(int(v_['1']),int(v_['NDIM'])+1):
                    v_['RR'+str(I)+','+str(J)+','+str(K)] = (v_['RG'+str(I)+','+str(J)]*v_['R'+
                         str(K)+','+str(J)])
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for J in range(int(v_['1']),int(v_['NBAR'])+1):
            [iv,ix_,_] = s2mpj_ii('XAREA'+str(J),ix_)
            self.xnames=arrset(self.xnames,iv,'XAREA'+str(J))
        for I in range(int(v_['1']),int(v_['NDIM'])+1):
            [iv,ix_,_] = s2mpj_ii('DISPL'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'DISPL'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for J in range(int(v_['1']),int(v_['NBAR'])+1):
            [ig,ig_,_] = s2mpj_ii('WEIGHT',ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['XAREA'+str(J)]])
            valA = np.append(valA,float(v_['W'+str(J)]))
        for K in range(int(v_['1']),int(v_['NDIM'])+1):
            [ig,ig_,_] = s2mpj_ii('EQUIL'+str(K),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'EQUIL'+str(K))
        for I in range(int(v_['1']),int(v_['NDIM'])+1):
            [ig,ig_,_] = s2mpj_ii('STREN',ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'STREN')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['DISPL'+str(I)]])
            valA = np.append(valA,float(v_['Q'+str(I)]))
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
        for K in range(int(v_['1']),int(v_['NDIM'])+1):
            self.gconst = arrset(self.gconst,ig_['EQUIL'+str(K)],float(v_['P'+str(K)]))
        self.gconst = arrset(self.gconst,ig_['STREN'],float(v_['ALPHA']))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        for J in range(int(v_['1']),int(v_['NBAR'])+1):
            self.xlower[ix_['XAREA'+str(J)]] = 1.0
        for I in range(int(v_['1']),int(v_['NDIM'])+1):
            self.xlower[ix_['DISPL'+str(I)]] = -float('Inf')
            self.xupper[ix_['DISPL'+str(I)]] = +float('Inf')
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'U')
        elftv = loaset(elftv,it,1,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['1']),int(v_['NDIM'])+1):
            for J in range(int(v_['1']),int(v_['NBAR'])+1):
                ename = 'UX'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'en2PR')
                ielftype = arrset(ielftype,ie,iet_["en2PR"])
                self.x0 = np.zeros((self.n,1))
                vname = 'DISPL'+str(I)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                posev = np.where(elftv[ielftype[ie]]=='U')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'XAREA'+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                posev = np.where(elftv[ielftype[ie]]=='X')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['NDIM'])+1):
            for J in range(int(v_['1']),int(v_['NBAR'])+1):
                for K in range(int(v_['1']),int(v_['NDIM'])+1):
                    ig = ig_['EQUIL'+str(K)]
                    posel = len(self.grelt[ig])
                    self.grelt = loaset(self.grelt,ig,posel,ie_['UX'+str(I)+','+str(J)])
                    nlc = np.union1d(nlc,np.array([ig]))
                    self.grelw  = (
                          loaset(self.grelw,ig,posel,float(v_['RR'+str(I)+','+str(J)+','+str(K)])))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Objective function value corresponding to the global minimizer above
        self.objlower = 1.2287408808
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
        self.pbclass   = "C-CLQR2-MN-11-4"
        self.x0        = np.zeros((self.n,1))
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def en2PR(self, nargout,*args):

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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

