from s2mpjlib import *
class  KOWOSB(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : KOWOSB
#    *********
# 
#    A problem arising in the analysis of kinetic data for an enzyme
#    reaction, known under the name of Kowalik and Osborne problem
#    in 4 variables.
# 
#    Source:  Problem 15 in
#    J.J. More', B.S. Garbow and K.E. Hillstrom,
#    "Testing Unconstrained Optimization Software",
#    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
# 
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "SUR2-MN-4-0"
# 
#    This function  is a nonlinear least squares with 11 groups.  Each
#    group has a linear and a nonlinear element.
# 
#    Number of groups
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'KOWOSB'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'KOWOSB'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['M'] = 11
        v_['1'] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2mpj_ii('X1',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X1')
        [iv,ix_,_] = s2mpj_ii('X2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X2')
        [iv,ix_,_] = s2mpj_ii('X3',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X3')
        [iv,ix_,_] = s2mpj_ii('X4',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X4')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['M'])+1):
            [ig,ig_,_] = s2mpj_ii('G'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        pbm.gconst = arrset(pbm.gconst,ig_['G1'],float(0.1957))
        pbm.gconst = arrset(pbm.gconst,ig_['G2'],float(0.1947))
        pbm.gconst = arrset(pbm.gconst,ig_['G3'],float(0.1735))
        pbm.gconst = arrset(pbm.gconst,ig_['G4'],float(0.1600))
        pbm.gconst = arrset(pbm.gconst,ig_['G5'],float(0.0844))
        pbm.gconst = arrset(pbm.gconst,ig_['G6'],float(0.0627))
        pbm.gconst = arrset(pbm.gconst,ig_['G7'],float(0.0456))
        pbm.gconst = arrset(pbm.gconst,ig_['G8'],float(0.0342))
        pbm.gconst = arrset(pbm.gconst,ig_['G9'],float(0.0323))
        pbm.gconst = arrset(pbm.gconst,ig_['G10'],float(0.0235))
        pbm.gconst = arrset(pbm.gconst,ig_['G11'],float(0.0246))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.x0[ix_['X1']] = float(0.25)
        pb.x0[ix_['X2']] = float(0.39)
        pb.x0[ix_['X3']] = float(0.415)
        pb.x0[ix_['X4']] = float(0.39)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eKWO', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftv = loaset(elftv,it,3,'V4')
        elftp = []
        elftp = loaset(elftp,it,0,'U')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['1']),int(v_['M'])+1):
            ename = 'E'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eKWO')
            ielftype = arrset(ielftype, ie, iet_["eKWO"])
            vname = 'X1'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X2'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X3'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X4'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V4')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        posep = find(elftp[ielftype[ie]],lambda x:x=='U')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(4.0))
        ename = 'E2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        posep = find(elftp[ielftype[ie]],lambda x:x=='U')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(2.0))
        ename = 'E3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        posep = find(elftp[ielftype[ie]],lambda x:x=='U')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        ename = 'E4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        posep = find(elftp[ielftype[ie]],lambda x:x=='U')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.5))
        ename = 'E5'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        posep = find(elftp[ielftype[ie]],lambda x:x=='U')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.25))
        ename = 'E6'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        posep = find(elftp[ielftype[ie]],lambda x:x=='U')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.167))
        ename = 'E7'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        posep = find(elftp[ielftype[ie]],lambda x:x=='U')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.125))
        ename = 'E8'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        posep = find(elftp[ielftype[ie]],lambda x:x=='U')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.1))
        ename = 'E9'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        posep = find(elftp[ielftype[ie]],lambda x:x=='U')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.0833))
        ename = 'E10'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        posep = find(elftp[ielftype[ie]],lambda x:x=='U')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.0714))
        ename = 'E11'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        posep = find(elftp[ielftype[ie]],lambda x:x=='U')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.0624))
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
        for I in range(int(v_['1']),int(v_['M'])+1):
            ig = ig_['G'+str(I)]
            pbm.grftype = arrset(pbm.grftype,ig,'gL2')
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        pb.objlower = 0.0
#    Solution
# LO SOLTN               0.00102734
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "SUR2-MN-4-0"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eKWO(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        USQ = pbm.elpar[iel_][0]*pbm.elpar[iel_][0]
        B1 = USQ+pbm.elpar[iel_][0]*EV_[1]
        B2 = USQ+pbm.elpar[iel_][0]*EV_[2]+EV_[3]
        B2SQ = B2*B2
        B2CB = B2*B2SQ
        UV1 = pbm.elpar[iel_][0]*EV_[0]
        UB1 = pbm.elpar[iel_][0]*B1
        T1 = B1/B2SQ
        T2 = 2.0/B2CB
        f_   = EV_[0]*B1/B2
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = B1/B2
            g_[1] = UV1/B2
            g_[2] = -UV1*T1
            g_[3] = -EV_[0]*T1
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,1] = pbm.elpar[iel_][0]/B2
                H_[1,0] = H_[0,1]
                H_[0,2] = -UB1/B2SQ
                H_[2,0] = H_[0,2]
                H_[0,3] = -T1
                H_[3,0] = H_[0,3]
                H_[1,2] = -UV1*pbm.elpar[iel_][0]/B2SQ
                H_[2,1] = H_[1,2]
                H_[1,3] = -UV1/B2SQ
                H_[3,1] = H_[1,3]
                H_[2,2] = T2*UV1*UB1
                H_[2,3] = T2*UV1*B1
                H_[3,2] = H_[2,3]
                H_[3,3] = T2*EV_[0]*B1
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
                H_ = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

