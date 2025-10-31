from s2mpjlib import *
class  SBRYBND(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : SBRYBND
#    *********
#    Broyden banded system of nonlinear equations, considered in the
#    least square sense.
#    NB: scaled version of BRYBND
# 
#    Source: problem 31 in
#    J.J. More', B.S. Garbow and K.E. Hillstrom,
#    "Testing Unconstrained Optimization Software",
#    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
# 
#    See also Buckley#73 (p. 41) and Toint#18
# 
#    SIF input: Ph. Toint and Nick Gould, Nov 1997.
# 
#    classification = "C-CSUR2-AN-V-0"
# 
#    N is the number of equations and variables (variable).
# 
#           Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER
# IE N                   50             $-PARAMETER
# IE N                   100            $-PARAMETER
# IE N                   500            $-PARAMETER
# IE N                   1000           $-PARAMETER     original value
# IE N                   5000           $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'SBRYBND'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(10);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
        v_['ONE'] = 1.0
        v_['KAPPA1'] = 2.0
        v_['KAPPA2'] = 5.0
        v_['KAPPA3'] = 1.0
        v_['LB'] = 5
        v_['UB'] = 1
        v_['RN'] = float(v_['N'])
        v_['RN-1'] = -1+v_['RN']
        v_['SCAL'] = 12.0
        v_['1'] = 1
        v_['MLB'] = -1*v_['LB']
        v_['MUB'] = -1*v_['UB']
        v_['LB+1'] = 1+v_['LB']
        v_['N-UB'] = v_['N']+v_['MUB']
        v_['N-UB-1'] = -1+v_['N-UB']
        v_['-KAPPA3'] = -1.0*v_['KAPPA3']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['I-1'] = -1+I
            v_['RI-1'] = float(v_['I-1'])
            v_['RAT'] = v_['RI-1']/v_['RN-1']
            v_['ARG'] = v_['RAT']*v_['SCAL']
            v_['SCALE'+str(I)] = np.exp(v_['ARG'])
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['1']),int(v_['LB'])+1):
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            v_['I+UB'] = I+v_['UB']
            for J in range(int(v_['1']),int(v_['I-1'])+1):
                v_['KAP'] = v_['-KAPPA3']*v_['SCALE'+str(J)]
                [ig,ig_,_] = s2mpj_ii('G'+str(I),ig_)
                gtype = arrset(gtype,ig,'<>')
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(J)]])
                valA = np.append(valA,float(v_['KAP']))
            v_['KAP'] = v_['KAPPA1']*v_['SCALE'+str(I)]
            [ig,ig_,_] = s2mpj_ii('G'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(I)]])
            valA = np.append(valA,float(v_['KAP']))
            for J in range(int(v_['I+1']),int(v_['I+UB'])+1):
                v_['KAP'] = v_['-KAPPA3']*v_['SCALE'+str(J)]
                [ig,ig_,_] = s2mpj_ii('G'+str(I),ig_)
                gtype = arrset(gtype,ig,'<>')
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(J)]])
                valA = np.append(valA,float(v_['KAP']))
        for I in range(int(v_['LB+1']),int(v_['N-UB-1'])+1):
            v_['I-LB'] = I+v_['MLB']
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            v_['I+UB'] = I+v_['UB']
            for J in range(int(v_['I-LB']),int(v_['I-1'])+1):
                v_['KAP'] = v_['-KAPPA3']*v_['SCALE'+str(J)]
                [ig,ig_,_] = s2mpj_ii('G'+str(I),ig_)
                gtype = arrset(gtype,ig,'<>')
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(J)]])
                valA = np.append(valA,float(v_['KAP']))
            v_['KAP'] = v_['KAPPA1']*v_['SCALE'+str(I)]
            [ig,ig_,_] = s2mpj_ii('G'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(I)]])
            valA = np.append(valA,float(v_['KAP']))
            for J in range(int(v_['I+1']),int(v_['I+UB'])+1):
                v_['KAP'] = v_['-KAPPA3']*v_['SCALE'+str(J)]
                [ig,ig_,_] = s2mpj_ii('G'+str(I),ig_)
                gtype = arrset(gtype,ig,'<>')
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(J)]])
                valA = np.append(valA,float(v_['KAP']))
        for I in range(int(v_['N-UB']),int(v_['N'])+1):
            v_['I-LB'] = I+v_['MLB']
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            for J in range(int(v_['I-LB']),int(v_['I-1'])+1):
                v_['KAP'] = v_['-KAPPA3']*v_['SCALE'+str(J)]
                [ig,ig_,_] = s2mpj_ii('G'+str(I),ig_)
                gtype = arrset(gtype,ig,'<>')
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(J)]])
                valA = np.append(valA,float(v_['KAP']))
            v_['KAP'] = v_['KAPPA1']*v_['SCALE'+str(I)]
            [ig,ig_,_] = s2mpj_ii('G'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(I)]])
            valA = np.append(valA,float(v_['KAP']))
            for J in range(int(v_['I+1']),int(v_['N'])+1):
                v_['KAP'] = v_['-KAPPA3']*v_['SCALE'+str(J)]
                [ig,ig_,_] = s2mpj_ii('G'+str(I),ig_)
                gtype = arrset(gtype,ig,'<>')
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(J)]])
                valA = np.append(valA,float(v_['KAP']))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['DIV'] = v_['ONE']/v_['SCALE'+str(I)]
            self.x0[ix_['X'+str(I)]] = float(v_['DIV'])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'V')
        elftp = []
        elftp = loaset(elftp,it,0,'P')
        [it,iet_,_] = s2mpj_ii( 'eCB', iet_)
        elftv = loaset(elftv,it,0,'V')
        elftp = loaset(elftp,it,0,'P')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['1']),int(v_['N'])+1):
            ename = 'E'+str(I)
            [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSQ')
            ielftype = arrset(ielftype,ie,iet_["eSQ"])
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='P')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['SCALE'+str(I)]))
            ename = 'Q'+str(I)
            [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eCB')
            ielftype = arrset(ielftype,ie,iet_["eCB"])
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='P')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['SCALE'+str(I)]))
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
        for I in range(int(v_['1']),int(v_['LB'])+1):
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            v_['I+UB'] = I+v_['UB']
            for J in range(int(v_['1']),int(v_['I-1'])+1):
                ig = ig_['G'+str(I)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(J)])
                self.grelw = loaset(self.grelw,ig,posel,float(v_['-KAPPA3']))
            ig = ig_['G'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['Q'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,float(v_['KAPPA2']))
            for J in range(int(v_['I+1']),int(v_['I+UB'])+1):
                ig = ig_['G'+str(I)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(J)])
                self.grelw = loaset(self.grelw,ig,posel,float(v_['-KAPPA3']))
        for I in range(int(v_['LB+1']),int(v_['N-UB-1'])+1):
            v_['I-LB'] = I+v_['MLB']
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            v_['I+UB'] = I+v_['UB']
            for J in range(int(v_['I-LB']),int(v_['I-1'])+1):
                ig = ig_['G'+str(I)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['Q'+str(J)])
                self.grelw = loaset(self.grelw,ig,posel,float(v_['-KAPPA3']))
            ig = ig_['G'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,float(v_['KAPPA2']))
            for J in range(int(v_['I+1']),int(v_['I+UB'])+1):
                ig = ig_['G'+str(I)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(J)])
                self.grelw = loaset(self.grelw,ig,posel,float(v_['-KAPPA3']))
        for I in range(int(v_['N-UB']),int(v_['N'])+1):
            v_['I-LB'] = I+v_['MLB']
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            for J in range(int(v_['I-LB']),int(v_['I-1'])+1):
                ig = ig_['G'+str(I)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(J)])
                self.grelw = loaset(self.grelw,ig,posel,float(v_['-KAPPA3']))
            ig = ig_['G'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['Q'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,float(v_['KAPPA2']))
            for J in range(int(v_['I+1']),int(v_['N'])+1):
                ig = ig_['G'+str(I)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(J)])
                self.grelw = loaset(self.grelw,ig,posel,float(v_['-KAPPA3']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               0.0
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass   = "C-CSUR2-AN-V-0"
        self.objderlvl = 2

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSQ(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        PP = self.elpar[iel_][0]*self.elpar[iel_][0]
        f_   = PP*EV_[0]*EV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = PP*(EV_[0]+EV_[0])
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0*PP
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eCB(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        PP = self.elpar[iel_][0]*self.elpar[iel_][0]*self.elpar[iel_][0]
        f_   = PP*EV_[0]*EV_[0]*EV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 3.0*PP*EV_[0]*EV_[0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 6.0*PP*EV_[0]
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

