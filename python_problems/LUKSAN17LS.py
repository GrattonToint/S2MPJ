from s2mpjlib import *
class  LUKSAN17LS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : LUKSAN17LS
#    *********
# 
#    Problem 17 (sparse trigonometric) in the paper
# 
#      L. Luksan
#      Hybrid methods in large sparse nonlinear least squares
#      J. Optimization Theory & Applications 89(3) 575-595 (1996)
# 
#    SIF input: Nick Gould, June 2017.
# 
#    least-squares version
# 
#    classification = "SUR2-AN-V-0"
# 
#   seed for dimensions
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'LUKSAN17LS'

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
        v_['S'] = 49
        v_['N'] = 2*v_['S']
        v_['N'] = 2+v_['N']
        v_['M'] = 4*v_['S']
        v_['1'] = 1
        v_['2'] = 2
        v_['3'] = 3
        v_['4'] = 4
        v_['ONE'] = 1.0
        v_['Y1'] = 30.6
        v_['Y2'] = 72.2
        v_['Y3'] = 124.4
        v_['Y4'] = 187.4
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
        for I in range(int(v_['1']),int(v_['M'])+1):
            [ig,ig_,_] = s2mpj_ii('E'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        v_['K'] = 1
        for J in range(int(v_['1']),int(v_['S'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['E'+str(int(v_['K']))],float(v_['Y1']))
            v_['K'] = 1+v_['K']
            pbm.gconst = arrset(pbm.gconst,ig_['E'+str(int(v_['K']))],float(v_['Y2']))
            v_['K'] = 1+v_['K']
            pbm.gconst = arrset(pbm.gconst,ig_['E'+str(int(v_['K']))],float(v_['Y3']))
            v_['K'] = 1+v_['K']
            pbm.gconst = arrset(pbm.gconst,ig_['E'+str(int(v_['K']))],float(v_['Y4']))
            v_['K'] = 1+v_['K']
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower = np.zeros((pb.n,1))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        for I in range(int(v_['1']),int(v_['N'])+1,int(v_['4'])):
            pb.x0[ix_['X'+str(I)]] = float(-0.8)
        for I in range(int(v_['2']),int(v_['N'])+1,int(v_['4'])):
            pb.x0[ix_['X'+str(I)]] = float(1.2)
        for I in range(int(v_['3']),int(v_['N'])+1,int(v_['4'])):
            pb.x0[ix_['X'+str(I)]] = float(-1.2)
        for I in range(int(v_['4']),int(v_['N'])+1,int(v_['4'])):
            pb.x0[ix_['X'+str(I)]] = float(0.8)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eACOSX', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftp = []
        elftp = loaset(elftp,it,0,'A')
        [it,iet_,_] = s2mpj_ii( 'eASINX', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftp = loaset(elftp,it,0,'A')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for Q in range(int(v_['1']),int(v_['4'])+1):
            v_['RQ'] = float(Q)
            v_['RQ2'] = v_['RQ']*v_['RQ']
            v_['K'] = 1
            v_['I'] = 0
            for J in range(int(v_['1']),int(v_['S'])+1):
                v_['I+Q'] = v_['I']+Q
                for L in range(int(v_['1']),int(v_['4'])+1):
                    v_['RL'] = float(L)
                    v_['RL2'] = v_['RL']*v_['RL']
                    v_['A'] = v_['RL']*v_['RQ2']
                    v_['A'] = -1.0*v_['A']
                    ename = 'S'+str(int(v_['K']))+','+str(Q)
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    pbm.elftype = arrset(pbm.elftype,ie,'eASINX')
                    ielftype = arrset(ielftype, ie, iet_["eASINX"])
                    ename = 'S'+str(int(v_['K']))+','+str(Q)
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    vname = 'X'+str(int(v_['I+Q']))
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='X')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    ename = 'S'+str(int(v_['K']))+','+str(Q)
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    posep = find(elftp[ielftype[ie]],lambda x:x=='A')
                    pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['A']))
                    v_['A'] = v_['RL2']*v_['RQ']
                    ename = 'C'+str(int(v_['K']))+','+str(Q)
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    pbm.elftype = arrset(pbm.elftype,ie,'eACOSX')
                    ielftype = arrset(ielftype, ie, iet_["eACOSX"])
                    ename = 'C'+str(int(v_['K']))+','+str(Q)
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    vname = 'X'+str(int(v_['I+Q']))
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='X')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    ename = 'C'+str(int(v_['K']))+','+str(Q)
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    posep = find(elftp[ielftype[ie]],lambda x:x=='A')
                    pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['A']))
                    v_['K'] = 1+v_['K']
                v_['I'] = 2+v_['I']
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
        for K in range(int(v_['1']),int(v_['M'])+1):
            for Q in range(int(v_['1']),int(v_['4'])+1):
                ig = ig_['E'+str(K)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['S'+str(K)+','+str(Q)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
                posel = posel+1
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C'+str(K)+','+str(Q)])
                pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN                0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "SUR2-AN-V-0"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eASINX(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        ASINX = pbm.elpar[iel_][0]*np.sin(EV_[0])
        f_   = ASINX
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = pbm.elpar[iel_][0]*np.cos(EV_[0])
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = -ASINX
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eACOSX(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        ACOSX = pbm.elpar[iel_][0]*np.cos(EV_[0])
        f_   = ACOSX
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -pbm.elpar[iel_][0]*np.sin(EV_[0])
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = -ACOSX
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

