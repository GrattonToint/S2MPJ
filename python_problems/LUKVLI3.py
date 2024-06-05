from s2mpjlib import *
class  LUKVLI3(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : LUKVLI3
#    *********
# 
#    Source: Problem 5.3, the chained Powell singular function with 
#    simplified trigonometric exponential constraints, 
#    due to L. Luksan and J. Vlcek,
#    "Sparse and partially separable test problems for 
#    unconstrained and equality constrained optimization",
#    Technical Report 767, Inst. Computer Science, Academy of Sciences
#    of the Czech Republic, 182 07 Prague, Czech Republic, 1999
# 
#    Equality constraints changed to inequalities
# 
#    SIF input: Nick Gould, April 2001
# 
#    classification = "OOR2-AY-V-V"
# 
#    some useful parameters, including N, the number of variables.
# 
#           Alternative values for the SIF file parameters:
# IE N                   100            $-PARAMETER
# IE N                   1000           $-PARAMETER
# IE N                   10000          $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'LUKVLI3'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'LUKVLI3'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(10);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
# IE N                   100000         $-PARAMETER
        v_['1'] = 1
        v_['2'] = 2
        v_['3'] = 3
        v_['4'] = 4
        v_['5'] = 5
        v_['6'] = 6
        v_['N/2'] = int(np.fix(v_['N']/v_['2']))
        v_['N/2-1'] = -1+v_['N/2']
        v_['N-1'] = -1+v_['N']
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
        for I in range(int(v_['1']),int(v_['N/2-1'])+1):
            v_['2I'] = 2*I
            v_['2I-1'] = -1+v_['2I']
            v_['2I+1'] = 1+v_['2I']
            v_['2I+2'] = 2+v_['2I']
            [ig,ig_,_] = s2mpj_ii('OA'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(int(v_['2I-1']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['2I']))]
            pbm.A[ig,iv] = float(10.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('OB'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(int(v_['2I+1']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['2I+2']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('OC'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(int(v_['2I']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['2I+1']))]
            pbm.A[ig,iv] = float(-2.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('OD'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(int(v_['2I-1']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['2I+2']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C1',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C1')
        iv = ix_['X'+str(int(v_['2']))]
        pbm.A[ig,iv] = float(2.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('C2',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'C2')
        iv = ix_['X'+str(int(v_['N-1']))]
        pbm.A[ig,iv] = float(4.0)+pbm.A[ig,iv]
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
        pbm.gconst = arrset(pbm.gconst,ig_['C1'],float(5.0))
        pbm.gconst = arrset(pbm.gconst,ig_['C2'],float(3.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        for I in range(int(v_['1']),int(v_['N'])+1,int(v_['4'])):
            pb.x0[ix_['X'+str(I)]] = float(3.0)
        for I in range(int(v_['2']),int(v_['N'])+1,int(v_['4'])):
            pb.x0[ix_['X'+str(I)]] = float(-1.0)
        for I in range(int(v_['3']),int(v_['N'])+1,int(v_['4'])):
            pb.x0[ix_['X'+str(I)]] = float(0.0)
        for I in range(int(v_['4']),int(v_['N'])+1,int(v_['4'])):
            pb.x0[ix_['X'+str(I)]] = float(1.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eCUBE', iet_)
        elftv = loaset(elftv,it,0,'V')
        [it,iet_,_] = s2mpj_ii( 'eXEXP', iet_)
        elftv = loaset(elftv,it,0,'X1')
        elftv = loaset(elftv,it,1,'X2')
        [it,iet_,_] = s2mpj_ii( 'eSINASINB', iet_)
        elftv = loaset(elftv,it,0,'X1')
        elftv = loaset(elftv,it,1,'X2')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        ename = 'C11'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCUBE')
        ielftype = arrset(ielftype, ie, iet_["eCUBE"])
        vname = 'X'+str(int(v_['1']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'C12'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSINASINB')
        ielftype = arrset(ielftype, ie, iet_["eSINASINB"])
        vname = 'X'+str(int(v_['1']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X'+str(int(v_['2']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'C21'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXEXP')
        ielftype = arrset(ielftype, ie, iet_["eXEXP"])
        vname = 'X'+str(int(v_['N-1']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X'+str(int(v_['N']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gL2P',igt_)
        [it,igt_,_] = s2mpj_ii('gL2P',igt_)
        grftp = []
        grftp = loaset(grftp,it,0,'P')
        [it,igt_,_] = s2mpj_ii('gL4P',igt_)
        [it,igt_,_] = s2mpj_ii('gL4P',igt_)
        grftp = loaset(grftp,it,0,'P')
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        pbm.grpar   = []
        for I in range(int(v_['1']),int(v_['N/2-1'])+1):
            ig = ig_['OA'+str(I)]
            pbm.grftype = arrset(pbm.grftype,ig,'gL2P')
            posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P')
            pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(1.0))
            ig = ig_['OB'+str(I)]
            pbm.grftype = arrset(pbm.grftype,ig,'gL2P')
            posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P')
            pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(5.0))
            ig = ig_['OC'+str(I)]
            pbm.grftype = arrset(pbm.grftype,ig,'gL4P')
            posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P')
            pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(1.0))
            ig = ig_['OD'+str(I)]
            pbm.grftype = arrset(pbm.grftype,ig,'gL4P')
            posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P')
            pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(10.0))
        ig = ig_['C1']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C11'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(3.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C12'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_['C2']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C21'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN               1.15775E+01
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.cupper[np.arange(pb.nle)] = np.zeros((pb.nle,1))
        pb.clower[np.arange(pb.nle+pb.neq,pb.m)] = np.zeros((pb.nge,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "OOR2-AY-V-V"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eXEXP(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        EXPX = np.exp(EV_[0]-EV_[1])
        XEXPX = EV_[0]*EXPX
        f_   = XEXPX
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EXPX+XEXPX
            g_[1] = -XEXPX
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 2.0*EXPX+XEXPX
                H_[0,1] = -EXPX-XEXPX
                H_[1,0] = H_[0,1]
                H_[1,1] = XEXPX
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eSINASINB(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        SINA = np.sin(EV_[0]-EV_[1])
        SINB = np.sin(EV_[0]+EV_[1])
        COSA = np.cos(EV_[0]-EV_[1])
        COSB = np.cos(EV_[0]+EV_[1])
        f_   = SINA*SINB
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = SINA*COSB+COSA*SINB
            g_[1] = SINA*COSB-COSA*SINB
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 2.0*(COSA*COSB-SINA*SINB)
                H_[1,1] = -2.0*(COSA*COSB+SINA*SINB)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eCUBE(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]**3
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 3.0*EV_[0]**2
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 6.0*EV_[0]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gL2P(pbm,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= pbm.grpar[igr_][0]*GVAR_*GVAR_
        if nargout>1:
            g_ = 2.0*pbm.grpar[igr_][0]*GVAR_
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 2.0*pbm.grpar[igr_][0]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def gL4P(pbm,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= pbm.grpar[igr_][0]*GVAR_**4
        if nargout>1:
            g_ = 4.0*pbm.grpar[igr_][0]*GVAR_**3
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 12.0*pbm.grpar[igr_][0]*GVAR_**2
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

