from s2xlib import *
class  PFIT1LS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : PFIT1LS
#    *********
# 
#    The problem is to fit a model containing a pole, given data
#    for values, first and second derivatives at two distinct points.
#    This is a least-squares version of problem PFIT1.
# 
#    The problem is not convex.
# 
#    SIF input: Ph. Toint, March 1994.
#               Lower bound on H added, Nov 2002.
# 
#    classification = "SBR2-AN-3-0"
# 
#    Problem data
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'PFIT1LS'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'PFIT1LS'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['CF'] = -8.0
        v_['CG'] = -18.66666666
        v_['CH'] = -23.11111111
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2x_ii('A',ix_)
        pb.xnames=arrset(pb.xnames,iv,'A')
        [iv,ix_,_] = s2x_ii('R',ix_)
        pb.xnames=arrset(pb.xnames,iv,'R')
        [iv,ix_,_] = s2x_ii('H',ix_)
        pb.xnames=arrset(pb.xnames,iv,'H')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2x_ii('EF',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2x_ii('EG',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2x_ii('EH',ig_)
        gtype = arrset(gtype,ig,'<>')
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        pbm.gconst = arrset(pbm.gconst,ig_['EF'],float(v_['CF']))
        pbm.gconst = arrset(pbm.gconst,ig_['EG'],float(v_['CG']))
        pbm.gconst = arrset(pbm.gconst,ig_['EH'],float(v_['CH']))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower[ix_['H']] = -0.5
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.x0[ix_['A']] = float(1.0)
        pb.x0[ix_['R']] = float(0.0)
        pb.x0[ix_['H']] = float(1.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eT1', iet_)
        elftv = loaset(elftv,it,0,'AA')
        elftv = loaset(elftv,it,1,'RR')
        elftv = loaset(elftv,it,2,'XX')
        [it,iet_,_] = s2x_ii( 'eT2', iet_)
        elftv = loaset(elftv,it,0,'AA')
        elftv = loaset(elftv,it,1,'RR')
        elftv = loaset(elftv,it,2,'XX')
        [it,iet_,_] = s2x_ii( 'eT3', iet_)
        elftv = loaset(elftv,it,0,'AA')
        elftv = loaset(elftv,it,1,'RR')
        elftv = loaset(elftv,it,2,'XX')
        [it,iet_,_] = s2x_ii( 'eT4', iet_)
        elftv = loaset(elftv,it,0,'AA')
        elftv = loaset(elftv,it,1,'RR')
        elftv = loaset(elftv,it,2,'XX')
        [it,iet_,_] = s2x_ii( 'eT5', iet_)
        elftv = loaset(elftv,it,0,'AA')
        elftv = loaset(elftv,it,1,'RR')
        elftv = loaset(elftv,it,2,'XX')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        ename = 'EA'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eT3')
        ielftype = arrset(ielftype, ie, iet_["eT3"])
        vname = 'A'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='AA')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'R'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RR')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'H'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EB'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eT2')
        ielftype = arrset(ielftype, ie, iet_["eT2"])
        vname = 'A'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='AA')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'R'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RR')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'H'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EC'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eT1')
        ielftype = arrset(ielftype, ie, iet_["eT1"])
        vname = 'A'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='AA')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'R'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RR')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'H'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'ED'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eT4')
        ielftype = arrset(ielftype, ie, iet_["eT4"])
        vname = 'A'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='AA')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'R'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RR')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'H'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EE'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eT5')
        ielftype = arrset(ielftype, ie, iet_["eT5"])
        vname = 'A'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='AA')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'R'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RR')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'H'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2x_ii('gL2',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for ig in range(0,ngrp):
            pbm.grftype = arrset(pbm.grftype,ig,'gL2')
        ig = ig_['EF']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EA'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-0.5))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EC'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['ED'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        ig = ig_['EG']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EA'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EB'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_['EH']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EE'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "SBR2-AN-3-0"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eT1(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[1]*EV_[2]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]*EV_[2]
            g_[1] = EV_[0]*EV_[2]
            g_[2] = EV_[0]*EV_[1]
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = EV_[2]
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1]
                H_[2,0] = H_[0,2]
                H_[1,2] = EV_[0]
                H_[2,1] = H_[1,2]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eT2(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        A1 = EV_[0]+1.0
        Y = 1.0+EV_[2]
        LOGY = np.log(Y)
        C = Y**(-A1)
        CC = C/Y
        CCC = CC/Y
        B = 1.0-C
        BA = LOGY*C
        BX = A1*CC
        BAA = -LOGY*LOGY*C
        BAX = -LOGY*BX+CC
        BXX = -A1*(A1+1.0)*CCC
        ARX = EV_[0]*EV_[1]*EV_[2]
        f_   = ARX*B
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]*EV_[2]*B+ARX*BA
            g_[1] = EV_[0]*EV_[2]*B
            g_[2] = EV_[0]*EV_[1]*B+ARX*BX
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = 2.0*EV_[1]*EV_[2]*BA+ARX*BAA
                H_[0,1] = EV_[2]*B+EV_[0]*EV_[2]*BA
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1]*B+EV_[1]*EV_[2]*BX+EV_[0]*EV_[1]*BA+ARX*BAX
                H_[2,0] = H_[0,2]
                H_[1,2] = EV_[0]*B+EV_[0]*EV_[2]*BX
                H_[2,1] = H_[1,2]
                H_[2,2] = 2.0*EV_[0]*EV_[1]*BX+ARX*BXX
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eT3(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*(EV_[0]+1.0)*EV_[1]*EV_[2]*EV_[2]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = (2.0*EV_[0]+1.0)*EV_[1]*EV_[2]*EV_[2]
            g_[1] = EV_[0]*(EV_[0]+1.0)*EV_[2]*EV_[2]
            g_[2] = 2.0*EV_[0]*(EV_[0]+1.0)*EV_[1]*EV_[2]
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = 2.0*EV_[1]*EV_[2]*EV_[2]
                H_[0,1] = (2.0*EV_[0]+1.0)*EV_[2]*EV_[2]
                H_[1,0] = H_[0,1]
                H_[0,2] = 2.0*(2.0*EV_[0]+1.0)*EV_[1]*EV_[2]
                H_[2,0] = H_[0,2]
                H_[1,2] = 2.0*EV_[0]*(EV_[0]+1.0)*EV_[2]
                H_[2,1] = H_[1,2]
                H_[2,2] = 2.0*EV_[0]*(EV_[0]+1.0)*EV_[1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eT4(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        Y = 1.0+EV_[2]
        LOGY = np.log(Y)
        C = Y**(-EV_[0])
        CC = C/Y
        CCC = CC/Y
        B = 1.0-C
        BA = LOGY*C
        BX = EV_[0]*CC
        BAA = -LOGY*LOGY*C
        BAX = -LOGY*BX+CC
        BXX = -EV_[0]*(EV_[0]+1.0)*CCC
        f_   = EV_[1]*B
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]*BA
            g_[1] = B
            g_[2] = EV_[1]*BX
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = EV_[1]*BAA
                H_[0,1] = BA
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1]*BAX
                H_[2,0] = H_[0,2]
                H_[1,2] = BX
                H_[2,1] = H_[1,2]
                H_[2,2] = EV_[1]*BXX
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eT5(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        A1 = EV_[0]+2.0
        Y = 1.0+EV_[2]
        LOGY = np.log(Y)
        C = Y**(-A1)
        CC = C/Y
        CCC = CC/Y
        B = 1.0-C
        BA = LOGY*C
        BX = A1*CC
        BAA = -LOGY*LOGY*C
        BAX = -LOGY*BX+CC
        BXX = -A1*(A1+1.0)*CCC
        D = EV_[0]*(EV_[0]+1.0)*EV_[1]*EV_[2]*EV_[2]
        DA = (2.0*EV_[0]+1.0)*EV_[1]*EV_[2]*EV_[2]
        DR = EV_[0]*(EV_[0]+1.0)*EV_[2]*EV_[2]
        DX = 2.0*EV_[0]*(EV_[0]+1.0)*EV_[1]*EV_[2]
        DAA = 2.0*EV_[1]*EV_[2]*EV_[2]
        DAR = (2.0*EV_[0]+1.0)*EV_[2]*EV_[2]
        DAX = 2.0*(2.0*EV_[0]+1.0)*EV_[1]*EV_[2]
        DRX = 2.0*EV_[0]*(EV_[0]+1.0)*EV_[2]
        DXX = 2.0*EV_[0]*(EV_[0]+1.0)*EV_[1]
        f_   = D*B
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = DA*B+D*BA
            g_[1] = DR*B
            g_[2] = DX*B+D*BX
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = DAA*B+2.0*DA*BA+D*BAA
                H_[0,1] = DAR*B+DR*BA
                H_[1,0] = H_[0,1]
                H_[0,2] = DAX*B+DA*BX+DX*BA+D*BAX
                H_[2,0] = H_[0,2]
                H_[1,2] = DRX*B+DR*BX
                H_[2,1] = H_[1,2]
                H_[2,2] = DXX*B+2.0*DX*BX+D*BXX
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

