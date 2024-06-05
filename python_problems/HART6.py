from s2mpjlib import *
class  HART6(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HART6
#    *********
# 
#    Source: Hartman problem 6 in
#    L. C. W. Dixon and G. P. Szego (Eds.)
#    Towards Global Optimization
#    North Holland, 1975.
#    Paper 9, page 163.
# 
#    SIF input: A.R. Conn May 1995
# 
#    classification = "OBR2-AN-6-0"
# 
#    Number of variables - constraints
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HART6'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'HART6'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 6
        v_['1'] = 1
        v_['ONE'] = 1
        v_['NN'] = 6
        v_['L'] = 4
        v_['C1'] = 1.0
        v_['C2'] = 1.2
        v_['C3'] = 3.0
        v_['C4'] = 3.2
        v_['A1,1'] = 10.0
        v_['A2,1'] = 0.05
        v_['A3,1'] = 3.0
        v_['A4,1'] = 17.0
        v_['A1,2'] = 0.05
        v_['A2,2'] = 10.0
        v_['A3,2'] = 3.5
        v_['A4,2'] = 8.0
        v_['A1,3'] = 17.0
        v_['A2,3'] = 17.0
        v_['A3,3'] = 1.7
        v_['A4,3'] = 0.05
        v_['A1,4'] = 3.5
        v_['A2,4'] = 0.1
        v_['A3,4'] = 10.0
        v_['A4,4'] = 10.0
        v_['A1,5'] = 1.7
        v_['A2,5'] = 8.0
        v_['A3,5'] = 17.0
        v_['A4,5'] = 0.1
        v_['A1,6'] = 8.0
        v_['A2,6'] = 14.0
        v_['A3,6'] = 8.0
        v_['A4,6'] = 14.0
        v_['P1,1'] = 0.1312
        v_['P2,1'] = 0.2329
        v_['P3,1'] = 0.2348
        v_['P4,1'] = 0.4047
        v_['P1,2'] = 0.1696
        v_['P2,2'] = 0.4135
        v_['P3,2'] = 0.1451
        v_['P4,2'] = 0.8828
        v_['P1,3'] = 0.5569
        v_['P2,3'] = 0.8307
        v_['P3,3'] = 0.3522
        v_['P4,3'] = 0.8732
        v_['P1,4'] = 0.0124
        v_['P2,4'] = 0.3736
        v_['P3,4'] = 0.2883
        v_['P4,4'] = 0.5743
        v_['P1,5'] = 0.8283
        v_['P2,5'] = 0.1004
        v_['P3,5'] = 0.3047
        v_['P4,5'] = 0.1091
        v_['P1,6'] = 0.5886
        v_['P2,6'] = 0.9991
        v_['P3,6'] = 0.6650
        v_['P4,6'] = 0.0381
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['L'])+1):
            [ig,ig_,_] = s2mpj_ii('OBJ'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            pbm.gscale = arrset(pbm.gscale,ig,float(-1.0))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),0.0)
        pb.xupper = np.full((pb.n,1),1.0)
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(0.2))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftp = []
        elftp = loaset(elftp,it,0,'PIJ')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['1']),int(v_['L'])+1):
            for J in range(int(v_['1']),int(v_['NN'])+1):
                ename = 'E'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
                ielftype = arrset(ielftype, ie, iet_["eSQ"])
                vname = 'X'+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,0.2)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                posep = find(elftp[ielftype[ie]],lambda x:x=='PIJ')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['P'+str(I)+','+str(J)]))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gNEXP',igt_)
        [it,igt_,_] = s2mpj_ii('gNEXP',igt_)
        grftp = []
        grftp = loaset(grftp,it,0,'CI')
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        pbm.grpar   = []
        for I in range(int(v_['1']),int(v_['L'])+1):
            ig = ig_['OBJ'+str(I)]
            pbm.grftype = arrset(pbm.grftype,ig,'gNEXP')
            for J in range(int(v_['1']),int(v_['NN'])+1):
                ig = ig_['OBJ'+str(I)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(I)+','+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['A'+str(I)+','+str(J)]))
            ig = ig_['OBJ'+str(I)]
            posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='CI')
            pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(v_['C'+str(I)]))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               -3.32288689158
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "OBR2-AN-6-0"
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
        f_   = (EV_[0]-pbm.elpar[iel_][0])*(EV_[0]-pbm.elpar[iel_][0])
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*(EV_[0]-pbm.elpar[iel_][0])
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gNEXP(pbm,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= pbm.grpar[igr_][0]*np.exp(-GVAR_)
        if nargout>1:
            g_ = -pbm.grpar[igr_][0]*np.exp(-GVAR_)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = pbm.grpar[igr_][0]*np.exp(-GVAR_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

