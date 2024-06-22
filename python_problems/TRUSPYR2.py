from s2mpjlib import *
class  TRUSPYR2(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    This is a structural optimization problem.
#    The problem is to minimize the weight of a given
#    8-bar truss structure formed as a pyramid for a given external load.
#    There are upper bounds on the normal stresses in the
#    bars and lower bounds on the cross-sectional areas of the bars.
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
#    classification = "LQR2-MN-11-11"
# 
#    Number of bars
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'TRUSPYR2'

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
        v_['NBAR'] = 8
        v_['NDIM'] = 3
        v_['1'] = 1
        v_['2'] = 2
        v_['NBAR/2'] = int(np.fix(v_['NBAR']/v_['2']))
        v_['8.0'] = 8.0
        v_['SQRT17'] = np.sqrt(17.0)
        v_['SQRT18'] = np.sqrt(18.0)
        v_['P1'] = 40.0
        v_['P2'] = 20.0
        v_['P3'] = 200.0
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
        pb.xnames = np.array([])
        pb.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for J in range(int(v_['1']),int(v_['NBAR'])+1):
            [iv,ix_,_] = s2mpj_ii('XAREA'+str(J),ix_)
            pb.xnames=arrset(pb.xnames,iv,'XAREA'+str(J))
        for I in range(int(v_['1']),int(v_['NDIM'])+1):
            [iv,ix_,_] = s2mpj_ii('DISPL'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'DISPL'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for J in range(int(v_['1']),int(v_['NBAR'])+1):
            [ig,ig_,_] = s2mpj_ii('WEIGHT',ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['XAREA'+str(J)]
            pbm.A[ig,iv] = float(v_['W'+str(J)])+pbm.A[ig,iv]
        for K in range(int(v_['1']),int(v_['NDIM'])+1):
            [ig,ig_,_] = s2mpj_ii('EQUIL'+str(K),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'EQUIL'+str(K))
        for I in range(int(v_['1']),int(v_['NDIM'])+1):
            for J in range(int(v_['1']),int(v_['NBAR'])+1):
                [ig,ig_,_] = s2mpj_ii('STRES'+str(J),ig_)
                gtype = arrset(gtype,ig,'<=')
                cnames = arrset(cnames,ig,'STRES'+str(J))
                iv = ix_['DISPL'+str(I)]
                pbm.A[ig,iv] = float(v_['R'+str(I)+','+str(J)])+pbm.A[ig,iv]
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
        for K in range(int(v_['1']),int(v_['NDIM'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['EQUIL'+str(K)],float(v_['P'+str(K)]))
        for J in range(int(v_['1']),int(v_['NBAR'])+1):
            pbm.gconst  = (
                  arrset(pbm.gconst,ig_['STRES'+str(J)],float(v_['STRUP'+str(J)])))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),float('inf'))
        for J in range(int(v_['1']),int(v_['NBAR'])+1):
            pb.xlower[ix_['XAREA'+str(J)]] = 1.0
        for I in range(int(v_['1']),int(v_['NDIM'])+1):
            pb.xlower[ix_['DISPL'+str(I)]] = -float('Inf')
            pb.xupper[ix_['DISPL'+str(I)]] = +float('Inf')
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'U')
        elftv = loaset(elftv,it,1,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['1']),int(v_['NDIM'])+1):
            for J in range(int(v_['1']),int(v_['NBAR'])+1):
                ename = 'UX'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
                ielftype = arrset(ielftype, ie, iet_["en2PR"])
                pb.x0 = np.zeros((pb.n,1))
                vname = 'DISPL'+str(I)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='U')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'XAREA'+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='X')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['NDIM'])+1):
            for J in range(int(v_['1']),int(v_['NBAR'])+1):
                for K in range(int(v_['1']),int(v_['NDIM'])+1):
                    ig = ig_['EQUIL'+str(K)]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['UX'+str(I)+','+str(J)])
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw  = (
                          loaset(pbm.grelw,ig,posel,float(v_['RR'+str(I)+','+str(J)+','+str(K)])))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Objective function value corresponding to the local minimizer above
        pb.objlower = 1.2287408808
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
        pb.pbclass = "LQR2-MN-11-11"
        pb.x0          = np.zeros((pb.n,1))
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def en2PR(pbm,nargout,*args):

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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

