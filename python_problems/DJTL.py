from s2mpjlib import *
class  DJTL(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DJTL
#    *********
# 
#    Source: modified version of problem 19 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
#    that is meant to simulate the Lagrangian barrier objective function
#    for particular values of the shifts and multipliers
# 
#    SIF input: A.R. Conn August 1993
# 
#    classification = "OUR2-AN-2-0"
# 
#    Define multipliers and shifts
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'DJTL'

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
        v_['LL1'] = 1.0
        v_['LL2'] = 1.0
        v_['LL3'] = 1.0
        v_['LL4'] = 1.0
        v_['SL1'] = 1.0
        v_['SL2'] = 1.0
        v_['SL3'] = 1.0
        v_['SL4'] = 1.0
        v_['LU1'] = 1.0
        v_['LU2'] = 1.0
        v_['LU3'] = 1.0
        v_['LU4'] = 1.0
        v_['SU1'] = 1.0
        v_['SU2'] = 1.0
        v_['SU3'] = 1.0
        v_['SU4'] = 1.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        pb.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2mpj_ii('X1',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X1')
        [iv,ix_,_] = s2mpj_ii('X2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X2')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('CONU1',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('CONL1',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('CONU2',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('CONL2',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('BNDU1',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('BNDL1',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('BNDU2',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X2']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('BNDL2',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        pbm.gconst = arrset(pbm.gconst,ig_['CONU1'],float(-200.0))
        pbm.gconst = arrset(pbm.gconst,ig_['CONL1'],float(100.0))
        pbm.gconst = arrset(pbm.gconst,ig_['CONU2'],float(0.0))
        pbm.gconst = arrset(pbm.gconst,ig_['CONL2'],float(-82.81))
        pbm.gconst = arrset(pbm.gconst,ig_['BNDU1'],float(-100.0))
        pbm.gconst = arrset(pbm.gconst,ig_['BNDL1'],float(13.0))
        pbm.gconst = arrset(pbm.gconst,ig_['BNDU2'],float(-100.0))
        pbm.gconst = arrset(pbm.gconst,ig_['BNDL2'],float(0.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower = np.zeros((pb.n,1))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.x0[ix_['X1']] = float(15.0)
        pb.x0[ix_['X2']] = float(6.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eCBm10', iet_)
        elftv = loaset(elftv,it,0,'V1')
        [it,iet_,_] = s2mpj_ii( 'eCBm20', iet_)
        elftv = loaset(elftv,it,0,'V1')
        [it,iet_,_] = s2mpj_ii( 'eSQm5', iet_)
        elftv = loaset(elftv,it,0,'V1')
        [it,iet_,_] = s2mpj_ii( 'eSQm6', iet_)
        elftv = loaset(elftv,it,0,'V1')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        ename = 'E1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCBm10')
        ielftype = arrset(ielftype, ie, iet_["eCBm10"])
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCBm20')
        ielftype = arrset(ielftype, ie, iet_["eCBm20"])
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSQm5')
        ielftype = arrset(ielftype, ie, iet_["eSQm5"])
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSQm5')
        ielftype = arrset(ielftype, ie, iet_["eSQm5"])
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E5'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSQm6')
        ielftype = arrset(ielftype, ie, iet_["eSQm6"])
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gLOG',igt_)
        [it,igt_,_] = s2mpj_ii('gLOG',igt_)
        grftp = []
        grftp = loaset(grftp,it,0,'P1')
        [it,igt_,_] = s2mpj_ii('gLOG',igt_)
        grftp = loaset(grftp,it,1,'P2')
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        pbm.grpar   = []
        ig = ig_['OBJ']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E1'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E2'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_['CONL1']
        pbm.grftype = arrset(pbm.grftype,ig,'gLOG')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E3'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E4'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P1')
        pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(v_['SL1']))
        posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P2')
        pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(v_['LL1']))
        ig = ig_['CONU1']
        pbm.grftype = arrset(pbm.grftype,ig,'gLOG')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E3'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E4'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P1')
        pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(v_['SU1']))
        posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P2')
        pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(v_['LU1']))
        ig = ig_['CONL2']
        pbm.grftype = arrset(pbm.grftype,ig,'gLOG')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E4'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E5'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P1')
        pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(v_['SL2']))
        posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P2')
        pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(v_['LL2']))
        ig = ig_['CONU2']
        pbm.grftype = arrset(pbm.grftype,ig,'gLOG')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E4'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E5'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P1')
        pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(v_['SU2']))
        posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P2')
        pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(v_['LU2']))
        ig = ig_['BNDL1']
        pbm.grftype = arrset(pbm.grftype,ig,'gLOG')
        posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P1')
        pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(v_['SL3']))
        posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P2')
        pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(v_['LL3']))
        ig = ig_['BNDU1']
        pbm.grftype = arrset(pbm.grftype,ig,'gLOG')
        posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P1')
        pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(v_['SU3']))
        posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P2')
        pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(v_['LU3']))
        ig = ig_['BNDL2']
        pbm.grftype = arrset(pbm.grftype,ig,'gLOG')
        posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P1')
        pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(v_['SL4']))
        posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P2')
        pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(v_['LL4']))
        ig = ig_['BNDU2']
        pbm.grftype = arrset(pbm.grftype,ig,'gLOG')
        posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P1')
        pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(v_['SU4']))
        posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P2')
        pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(v_['LU4']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               -8951.54472
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "OUR2-AN-2-0"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eCBm10(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        DIF = EV_[0]-10.0
        f_   = DIF**3
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 3.0*DIF*DIF
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 6.0*DIF
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eCBm20(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        DIF = EV_[0]-20.0
        f_   = DIF**3
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 3.0*DIF*DIF
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 6.0*DIF
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eSQm5(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        DIF = EV_[0]-5.0
        f_   = DIF**2
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*DIF
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
    def eSQm6(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        DIF = EV_[0]-6.0
        f_   = DIF**2
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*DIF
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
    def gLOG(pbm,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        APP1 = GVAR_+pbm.grpar[igr_][0]
        P1P2 = pbm.grpar[igr_][0]*pbm.grpar[igr_][1]
        ARG0 = APP1<=0.0
        BIG = 1.0000e+10
        if ARG0!=0:
            FF = BIG*GVAR_**2
        if ARG0==0:
            FF = -P1P2*np.log(APP1)
        if ARG0!=0:
            GG = 2.0*BIG*GVAR_
        if ARG0==0:
            GG = -P1P2/APP1
        if ARG0!=0:
            HH = 2.0*BIG
        if ARG0==0:
            HH = P1P2/APP1**2
        f_= FF
        if nargout>1:
            g_ = GG
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = HH
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

