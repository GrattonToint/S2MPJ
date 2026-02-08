from s2mpjlib import *
class  VAREIGVL(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : VAREIGVL
#    *********
# 
#    The variational eigenvalue by Auchmuty.
#    This problems features a banded matrix of bandwidth 2M+1 = 9.
# 
#    This problem has N least-squares groups, each having a linear part
#    only and N nonlinear elements,
#    plus a least q-th power group having N nonlinear elements.
# 
#    Source: problem 1 in
#    J.J. More',
#    "A collection of nonlinear model problems"
#    Proceedings of the AMS-SIAM Summer seminar on the Computational
#    Solution of Nonlinear Systems of Equations, Colorado, 1988.
#    Argonne National Laboratory MCS-P60-0289, 1989.
# 
#    SIF input: Ph. Toint, Dec 1989.
#               correction by Ph. Shott, January, 1995
#               and Nick Gould, December, 2019, May 2024
# 
#    classification = "C-COUR2-AN-V-0"
# 
#    Number of variables -1 (variable)
# 
#           Alternative values for the SIF file parameters:
# IE N                   19             $-PARAMETER
# IE N                   49             $-PARAMETER     original value
# IE N                   99             $-PARAMETER
# IE N                   499            $-PARAMETER
# IE N                   999            $-PARAMETER
# IE N                   4999           $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'VAREIGVL'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(19);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
# IE M                   4              $-PARAMETER  .le. N
# IE M                   5              $-PARAMETER  .le. N
        if nargin<2:
            v_['M'] = int(6);  #  SIF file default value
        else:
            v_['M'] = int(args[1])
        if nargin<3:
            v_['Q'] = float(1.5);  #  SIF file default value
        else:
            v_['Q'] = float(args[2])
        v_['1'] = 1
        v_['-1.0'] = -1.0
        v_['N+1'] = 1+v_['N']
        v_['-M'] = -1*v_['M']
        v_['M+1'] = 1+v_['M']
        v_['N-M'] = v_['N']+v_['-M']
        v_['N-M+1'] = 1+v_['N-M']
        v_['N2'] = v_['N']*v_['N']
        v_['RN2'] = float(v_['N2'])
        v_['-1/N2'] = v_['-1.0']/v_['RN2']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I))
        [iv,ix_,_] = s2mpj_ii('MU',ix_)
        self.xnames=arrset(self.xnames,iv,'MU')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['1']),int(v_['M'])+1):
            v_['RI'] = float(I)
            v_['-I'] = -1.0*v_['RI']
            v_['I+M'] = I+v_['M']
            for J in range(int(v_['1']),int(v_['I+M'])+1):
                v_['RJ'] = float(J)
                v_['IJ'] = v_['RI']*v_['RJ']
                v_['SIJ'] = np.sin(v_['IJ'])
                v_['J-I'] = v_['RJ']+v_['-I']
                v_['J-ISQ'] = v_['J-I']*v_['J-I']
                v_['ARG'] = v_['J-ISQ']*v_['-1/N2']
                v_['EX'] = np.exp(v_['ARG'])
                v_['AIJ'] = v_['SIJ']*v_['EX']
                [ig,ig_,_] = s2mpj_ii('G'+str(I),ig_)
                gtype = arrset(gtype,ig,'<>')
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(J)]])
                valA = np.append(valA,float(v_['AIJ']))
        for I in range(int(v_['M+1']),int(v_['N-M'])+1):
            v_['RI'] = float(I)
            v_['-I'] = -1.0*v_['RI']
            v_['I-M'] = I+v_['-M']
            v_['I+M'] = I+v_['M']
            for J in range(int(v_['I-M']),int(v_['I+M'])+1):
                v_['RJ'] = float(J)
                v_['IJ'] = v_['RI']*v_['RJ']
                v_['SIJ'] = np.sin(v_['IJ'])
                v_['J-I'] = v_['RJ']+v_['-I']
                v_['J-ISQ'] = v_['J-I']*v_['J-I']
                v_['ARG'] = v_['J-ISQ']*v_['-1/N2']
                v_['EX'] = np.exp(v_['ARG'])
                v_['AIJ'] = v_['SIJ']*v_['EX']
                [ig,ig_,_] = s2mpj_ii('G'+str(I),ig_)
                gtype = arrset(gtype,ig,'<>')
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(J)]])
                valA = np.append(valA,float(v_['AIJ']))
        for I in range(int(v_['N-M+1']),int(v_['N'])+1):
            v_['RI'] = float(I)
            v_['-I'] = -1.0*v_['RI']
            v_['I-M'] = I+v_['-M']
            for J in range(int(v_['I-M']),int(v_['N'])+1):
                v_['RJ'] = float(J)
                v_['IJ'] = v_['RI']*v_['RJ']
                v_['SIJ'] = np.sin(v_['IJ'])
                v_['J-I'] = v_['RJ']+v_['-I']
                v_['J-ISQ'] = v_['J-I']*v_['J-I']
                v_['ARG'] = v_['J-ISQ']*v_['-1/N2']
                v_['EX'] = np.exp(v_['ARG'])
                v_['AIJ'] = v_['SIJ']*v_['EX']
                [ig,ig_,_] = s2mpj_ii('G'+str(I),ig_)
                gtype = arrset(gtype,ig,'<>')
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(J)]])
                valA = np.append(valA,float(v_['AIJ']))
        [ig,ig_,_] = s2mpj_ii('G'+str(int(v_['N+1'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        selfnob      = ngrp
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(1.0))
        self.x0[ix_['MU']] = float(0.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'M')
        elftv = loaset(elftv,it,1,'X')
        [it,iet_,_] = s2mpj_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['1']),int(v_['N'])+1):
            ename = 'P'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_["en2PR"])
            vname = 'MU'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(1.0))
            posev = np.where(elftv[ielftype[ie]]=='M')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(1.0))
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'S'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSQ')
            ielftype = arrset(ielftype,ie,iet_["eSQ"])
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(1.0))
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gLQ',igt_)
        [it,igt_,_] = s2mpj_ii('gLQ',igt_)
        grftp = []
        grftp = loaset(grftp,it,0,'POWER')
        [it,igt_,_] = s2mpj_ii('gLQ2',igt_)
        [it,igt_,_] = s2mpj_ii('gLQ2',igt_)
        grftp = loaset(grftp,it,0,'POWER')
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        self.grpar   = []
        for I in range(int(v_['1']),int(v_['N'])+1):
            ig = ig_['G'+str(I)]
            self.grftype = arrset(self.grftype,ig,'gLQ')
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['P'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
            posgp = np.where(grftp[igt_[self.grftype[ig]]]=='POWER')[0]
            self.grpar =loaset(self.grpar,ig,posgp[0],float(2.0))
        ig = ig_['G'+str(int(v_['N+1']))]
        self.grftype = arrset(self.grftype,ig,'gLQ2')
        for I in range(int(v_['1']),int(v_['N'])+1):
            ig = ig_['G'+str(int(v_['N+1']))]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['S'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['G'+str(int(v_['N+1']))]
        posgp = np.where(grftp[igt_[self.grftype[ig]]]=='POWER')[0]
        self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['Q']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 0.0
#    Solution
# LO SOLTN               0.0
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass   = "C-COUR2-AN-V-0"
        self.objderlvl = 2

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

    @staticmethod
    def eSQ(self, nargout,*args):

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

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gLQ(self,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        IPOWER = self.grpar[igr_][0]
        PM1 = IPOWER-1
        f_= GVAR_**IPOWER/self.grpar[igr_][0]
        if nargout>1:
            g_ = GVAR_**PM1
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = PM1*GVAR_**(IPOWER-2)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def gLQ2(self,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= GVAR_**self.grpar[igr_][0]/self.grpar[igr_][0]
        if nargout>1:
            g_ = GVAR_**(self.grpar[igr_][0]-1.0e0)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = (self.grpar[igr_][0]-1.0e0)*GVAR_**(self.grpar[igr_][0]-2.0e0)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

