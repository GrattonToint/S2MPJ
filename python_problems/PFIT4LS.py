from s2mpjlib import *
class  PFIT4LS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : PFIT4LS
#    *********
# 
#    The problem is to fit a model containing a pole, given data
#    for values, first and second derivatives at two distinct points.
#    This is a least-squares version of problem PFIT4.
# 
#    The problem is not convex.
# 
#    SIF input: Ph. Toint, March 1994.
#               Lower bound on H added, Nov 2002.
# 
#    classification = "C-CSBR2-AN-3-0"
# 
#    Problem data
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'PFIT4LS'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['CF'] = -98.96296296
        v_['CG'] = -216.0987654
        v_['CH'] = -239.6707818
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [iv,ix_,_] = s2mpj_ii('A',ix_)
        self.xnames=arrset(self.xnames,iv,'A')
        [iv,ix_,_] = s2mpj_ii('R',ix_)
        self.xnames=arrset(self.xnames,iv,'R')
        [iv,ix_,_] = s2mpj_ii('H',ix_)
        self.xnames=arrset(self.xnames,iv,'H')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('EF',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('EG',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('EH',ig_)
        gtype = arrset(gtype,ig,'<>')
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        selfnob      = ngrp
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        self.gconst = arrset(self.gconst,ig_['EF'],float(v_['CF']))
        self.gconst = arrset(self.gconst,ig_['EG'],float(v_['CG']))
        self.gconst = arrset(self.gconst,ig_['EH'],float(v_['CH']))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        self.xlower[ix_['H']] = -0.5
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.x0[ix_['A']] = float(1.0)
        self.x0[ix_['R']] = float(0.0)
        self.x0[ix_['H']] = float(1.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eT1', iet_)
        elftv = loaset(elftv,it,0,'AA')
        elftv = loaset(elftv,it,1,'RR')
        elftv = loaset(elftv,it,2,'XX')
        [it,iet_,_] = s2mpj_ii( 'eT2', iet_)
        elftv = loaset(elftv,it,0,'AA')
        elftv = loaset(elftv,it,1,'RR')
        elftv = loaset(elftv,it,2,'XX')
        [it,iet_,_] = s2mpj_ii( 'eT3', iet_)
        elftv = loaset(elftv,it,0,'AA')
        elftv = loaset(elftv,it,1,'RR')
        elftv = loaset(elftv,it,2,'XX')
        [it,iet_,_] = s2mpj_ii( 'eT4', iet_)
        elftv = loaset(elftv,it,0,'AA')
        elftv = loaset(elftv,it,1,'RR')
        elftv = loaset(elftv,it,2,'XX')
        [it,iet_,_] = s2mpj_ii( 'eT5', iet_)
        elftv = loaset(elftv,it,0,'AA')
        elftv = loaset(elftv,it,1,'RR')
        elftv = loaset(elftv,it,2,'XX')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        ename = 'EA'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eT3')
        ielftype = arrset(ielftype,ie,iet_["eT3"])
        vname = 'A'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='AA')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'R'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='RR')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'H'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='XX')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'EB'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eT2')
        ielftype = arrset(ielftype,ie,iet_["eT2"])
        vname = 'A'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='AA')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'R'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='RR')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'H'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='XX')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'EC'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eT1')
        ielftype = arrset(ielftype,ie,iet_["eT1"])
        vname = 'A'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='AA')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'R'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='RR')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'H'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='XX')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'ED'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eT4')
        ielftype = arrset(ielftype,ie,iet_["eT4"])
        vname = 'A'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='AA')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'R'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='RR')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'H'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='XX')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'EE'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eT5')
        ielftype = arrset(ielftype,ie,iet_["eT5"])
        vname = 'A'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='AA')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'R'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='RR')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'H'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='XX')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gL2',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for ig in range(0,ngrp):
            self.grftype = arrset(self.grftype,ig,'gL2')
        ig = ig_['EF']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['EA'])
        self.grelw = loaset(self.grelw,ig,posel,float(-0.5))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['EC'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['ED'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        ig = ig_['EG']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['EA'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['EB'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        ig = ig_['EH']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['EE'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 0.0
#    Solution at ( 1.0, 3.0 , 2.0 )
# LO SOLTN               0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass   = "C-CSBR2-AN-3-0"
        self.objderlvl = 2

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eT1(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0,0]*EV_[1,0]*EV_[2,0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1,0]*EV_[2,0]
            g_[1] = EV_[0,0]*EV_[2,0]
            g_[2] = EV_[0,0]*EV_[1,0]
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = EV_[2,0]
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1,0]
                H_[2,0] = H_[0,2]
                H_[1,2] = EV_[0,0]
                H_[2,1] = H_[1,2]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eT2(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        A1 = EV_[0,0]+1.0
        Y = 1.0+EV_[2,0]
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
        ARX = EV_[0,0]*EV_[1,0]*EV_[2,0]
        f_   = ARX*B
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1,0]*EV_[2,0]*B+ARX*BA
            g_[1] = EV_[0,0]*EV_[2,0]*B
            g_[2] = EV_[0,0]*EV_[1,0]*B+ARX*BX
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = 2.0*EV_[1,0]*EV_[2,0]*BA+ARX*BAA
                H_[0,1] = EV_[2,0]*B+EV_[0,0]*EV_[2,0]*BA
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1,0]*B+EV_[1,0]*EV_[2,0]*BX+EV_[0,0]*EV_[1,0]*BA+ARX*BAX
                H_[2,0] = H_[0,2]
                H_[1,2] = EV_[0,0]*B+EV_[0,0]*EV_[2,0]*BX
                H_[2,1] = H_[1,2]
                H_[2,2] = 2.0*EV_[0,0]*EV_[1,0]*BX+ARX*BXX
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eT3(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0,0]*(EV_[0,0]+1.0)*EV_[1,0]*EV_[2,0]*EV_[2,0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = (2.0*EV_[0,0]+1.0)*EV_[1,0]*EV_[2,0]*EV_[2,0]
            g_[1] = EV_[0,0]*(EV_[0,0]+1.0)*EV_[2,0]*EV_[2,0]
            g_[2] = 2.0*EV_[0,0]*(EV_[0,0]+1.0)*EV_[1,0]*EV_[2,0]
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = 2.0*EV_[1,0]*EV_[2,0]*EV_[2,0]
                H_[0,1] = (2.0*EV_[0,0]+1.0)*EV_[2,0]*EV_[2,0]
                H_[1,0] = H_[0,1]
                H_[0,2] = 2.0*(2.0*EV_[0,0]+1.0)*EV_[1,0]*EV_[2,0]
                H_[2,0] = H_[0,2]
                H_[1,2] = 2.0*EV_[0,0]*(EV_[0,0]+1.0)*EV_[2,0]
                H_[2,1] = H_[1,2]
                H_[2,2] = 2.0*EV_[0,0]*(EV_[0,0]+1.0)*EV_[1,0]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eT4(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        Y = 1.0+EV_[2,0]
        LOGY = np.log(Y)
        C = Y**(-EV_[0,0])
        CC = C/Y
        CCC = CC/Y
        B = 1.0-C
        BA = LOGY*C
        BX = EV_[0,0]*CC
        BAA = -LOGY*LOGY*C
        BAX = -LOGY*BX+CC
        BXX = -EV_[0,0]*(EV_[0,0]+1.0)*CCC
        f_   = EV_[1,0]*B
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1,0]*BA
            g_[1] = B
            g_[2] = EV_[1,0]*BX
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = EV_[1,0]*BAA
                H_[0,1] = BA
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1,0]*BAX
                H_[2,0] = H_[0,2]
                H_[1,2] = BX
                H_[2,1] = H_[1,2]
                H_[2,2] = EV_[1,0]*BXX
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eT5(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        A1 = EV_[0,0]+2.0
        Y = 1.0+EV_[2,0]
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
        D = EV_[0,0]*(EV_[0,0]+1.0)*EV_[1,0]*EV_[2,0]*EV_[2,0]
        DA = (2.0*EV_[0,0]+1.0)*EV_[1,0]*EV_[2,0]*EV_[2,0]
        DR = EV_[0,0]*(EV_[0,0]+1.0)*EV_[2,0]*EV_[2,0]
        DX = 2.0*EV_[0,0]*(EV_[0,0]+1.0)*EV_[1,0]*EV_[2,0]
        DAA = 2.0*EV_[1,0]*EV_[2,0]*EV_[2,0]
        DAR = (2.0*EV_[0,0]+1.0)*EV_[2,0]*EV_[2,0]
        DAX = 2.0*(2.0*EV_[0,0]+1.0)*EV_[1,0]*EV_[2,0]
        DRX = 2.0*EV_[0,0]*(EV_[0,0]+1.0)*EV_[2,0]
        DXX = 2.0*EV_[0,0]*(EV_[0,0]+1.0)*EV_[1,0]
        f_   = D*B
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
    def gL2(self,nargout,*args):

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

