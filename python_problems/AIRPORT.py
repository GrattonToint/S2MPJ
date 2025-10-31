from s2mpjlib import *
class  AIRPORT(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    This problem is concerned with the localisation of airports in Brazil.
#    We consider  m  balls in the real plane, whose centers are the coordinates
#    of some Brazilian  cities and whose  radius were chosen such that the balls are
#    disjoint. The problem is to find one point  (xi, yi) on  each ball, i=1,..,m,
#    such that  SUM(||(xi,yi) - (xj,yj)||)  is  minimum, where the sum involves all
#    the pairs (i,j) such that 1 <= i <= m, 1 <= j <= m and i <> j.
# 
#    For this problem instance, we have m =  42 cities and n = 84 points, 
#    i.e, 42 nonlinear inequalities constraints and 84 variables.
# 
#    Source:
#    Contribution from a LANCELOT user.
# 
#    SIF input : Rodrigo de Barros Nabholz & Maria Aparecida Diniz Ehrhardt
#                November 1994, DMA - IMECC- UNICAMP
#    Adaptation for CUTE: Ph. Toint, November 1994.
# 
#    classification = "C-CSQR2-MN-84-42"
# 
#    Problem data
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'AIRPORT'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 42
        v_['N-1'] = 41
        v_['1'] = 1
        v_['R1'] = 0.09
        v_['R2'] = 0.3
        v_['R3'] = 0.09
        v_['R4'] = 0.45
        v_['R5'] = 0.5
        v_['R6'] = 0.04
        v_['R7'] = 0.1
        v_['R8'] = 0.02
        v_['R9'] = 0.02
        v_['R10'] = 0.07
        v_['R11'] = 0.4
        v_['R12'] = 0.045
        v_['R13'] = 0.05
        v_['R14'] = 0.056
        v_['R15'] = 0.36
        v_['R16'] = 0.08
        v_['R17'] = 0.07
        v_['R18'] = 0.36
        v_['R19'] = 0.67
        v_['R20'] = 0.38
        v_['R21'] = 0.37
        v_['R22'] = 0.05
        v_['R23'] = 0.4
        v_['R24'] = 0.66
        v_['R25'] = 0.05
        v_['R26'] = 0.07
        v_['R27'] = 0.08
        v_['R28'] = 0.3
        v_['R29'] = 0.31
        v_['R30'] = 0.49
        v_['R31'] = 0.09
        v_['R32'] = 0.46
        v_['R33'] = 0.12
        v_['R34'] = 0.07
        v_['R35'] = 0.07
        v_['R36'] = 0.09
        v_['R37'] = 0.05
        v_['R38'] = 0.13
        v_['R39'] = 0.16
        v_['R40'] = 0.46
        v_['R41'] = 0.25
        v_['R42'] = 0.1
        v_['CX1'] = -6.3
        v_['CX2'] = -7.8
        v_['CX3'] = -9.0
        v_['CX4'] = -7.2
        v_['CX5'] = -5.7
        v_['CX6'] = -1.9
        v_['CX7'] = -3.5
        v_['CX8'] = -0.5
        v_['CX9'] = 1.4
        v_['CX10'] = 4.0
        v_['CX11'] = 2.1
        v_['CX12'] = 5.5
        v_['CX13'] = 5.7
        v_['CX14'] = 5.7
        v_['CX15'] = 3.8
        v_['CX16'] = 5.3
        v_['CX17'] = 4.7
        v_['CX18'] = 3.3
        v_['CX19'] = 0.0
        v_['CX20'] = -1.0
        v_['CX21'] = -0.4
        v_['CX22'] = 4.2
        v_['CX23'] = 3.2
        v_['CX24'] = 1.7
        v_['CX25'] = 3.3
        v_['CX26'] = 2.0
        v_['CX27'] = 0.7
        v_['CX28'] = 0.1
        v_['CX29'] = -0.1
        v_['CX30'] = -3.5
        v_['CX31'] = -4.0
        v_['CX32'] = -2.7
        v_['CX33'] = -0.5
        v_['CX34'] = -2.9
        v_['CX35'] = -1.2
        v_['CX36'] = -0.4
        v_['CX37'] = -0.1
        v_['CX38'] = -1.0
        v_['CX39'] = -1.7
        v_['CX40'] = -2.1
        v_['CX41'] = -1.8
        v_['CX42'] = 0.0
        v_['CY1'] = 8.0
        v_['CY2'] = 5.1
        v_['CY3'] = 2.0
        v_['CY4'] = 2.6
        v_['CY5'] = 5.5
        v_['CY6'] = 7.1
        v_['CY7'] = 5.9
        v_['CY8'] = 6.6
        v_['CY9'] = 6.1
        v_['CY10'] = 5.6
        v_['CY11'] = 4.9
        v_['CY12'] = 4.7
        v_['CY13'] = 4.3
        v_['CY14'] = 3.6
        v_['CY15'] = 4.1
        v_['CY16'] = 3.0
        v_['CY17'] = 2.4
        v_['CY18'] = 3.0
        v_['CY19'] = 4.7
        v_['CY20'] = 3.4
        v_['CY21'] = 2.3
        v_['CY22'] = 1.5
        v_['CY23'] = 0.5
        v_['CY24'] = -1.7
        v_['CY25'] = -2.0
        v_['CY26'] = -3.1
        v_['CY27'] = -3.5
        v_['CY28'] = -2.4
        v_['CY29'] = -1.3
        v_['CY30'] = 0.0
        v_['CY31'] = -1.7
        v_['CY32'] = -2.1
        v_['CY33'] = -0.4
        v_['CY34'] = -2.9
        v_['CY35'] = -3.4
        v_['CY36'] = -4.3
        v_['CY37'] = -5.2
        v_['CY38'] = -6.5
        v_['CY39'] = -7.5
        v_['CY40'] = -6.4
        v_['CY41'] = -5.1
        v_['CY42'] = 0.0
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
            [iv,ix_,_] = s2mpj_ii('Y'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'Y'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            v_['I+1'] = I+v_['1']
            for J in range(int(v_['I+1']),int(v_['N'])+1):
                [ig,ig_,_] = s2mpj_ii('OBJ1'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(I)]])
                valA = np.append(valA,float(1.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(J)]])
                valA = np.append(valA,float(-1.0))
                [ig,ig_,_] = s2mpj_ii('OBJ2'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['Y'+str(I)]])
                valA = np.append(valA,float(1.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['Y'+str(J)]])
                valA = np.append(valA,float(-1.0))
        for I in range(int(v_['1']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('CONS'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'CONS'+str(I))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        legrps = np.where(gtype=='<=')[0]
        eqgrps = np.where(gtype=='==')[0]
        gegrps = np.where(gtype=='>=')[0]
        self.nle = len(legrps)
        self.neq = len(eqgrps)
        self.nge = len(gegrps)
        self.m   = self.nle+self.neq+self.nge
        self.congrps = np.concatenate((legrps,eqgrps,gegrps))
        self.cnames = cnames[self.congrps]
        self.nob = ngrp-self.m
        self.objgrps = np.where(gtype=='<>')[0]
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        for I in range(int(v_['1']),int(v_['N'])+1):
            self.gconst = arrset(self.gconst,ig_['CONS'+str(I)],float(v_['R'+str(I)]))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        for I in range(int(v_['1']),int(v_['N'])+1):
            self.xlower[ix_['X'+str(I)]] = -10
            self.xupper[ix_['X'+str(I)]] = 10
            self.xlower[ix_['Y'+str(I)]] = -10
            self.xupper[ix_['Y'+str(I)]] = 10
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eDIFSQR', iet_)
        elftv = loaset(elftv,it,0,'V')
        elftp = []
        elftp = loaset(elftp,it,0,'W')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['1']),int(v_['N'])+1):
            ename = 'A'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eDIFSQR')
            ielftype = arrset(ielftype,ie,iet_["eDIFSQR"])
            self.x0 = np.zeros((self.n,1))
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='W')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['CX'+str(I)]))
        for I in range(int(v_['1']),int(v_['N'])+1):
            ename = 'B'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eDIFSQR')
            ielftype = arrset(ielftype,ie,iet_["eDIFSQR"])
            vname = 'Y'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='W')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['CY'+str(I)]))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gSQUARE',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            v_['I+1'] = I+v_['1']
            for J in range(int(v_['I+1']),int(v_['N'])+1):
                ig = ig_['OBJ1'+str(I)+','+str(J)]
                self.grftype = arrset(self.grftype,ig,'gSQUARE')
                ig = ig_['OBJ2'+str(I)+','+str(J)]
                self.grftype = arrset(self.grftype,ig,'gSQUARE')
        for I in range(int(v_['1']),int(v_['N'])+1):
            ig = ig_['CONS'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['A'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(1.0))
            posel = posel+1
            self.grelt = loaset(self.grelt,ig,posel,ie_['B'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = .0
#    Solution
# LO SOLTN              47952.695811
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.cupper[np.arange(self.nle)] = np.zeros((self.nle,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CSQR2-MN-84-42"
        self.x0        = np.zeros((self.n,1))
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eDIFSQR(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        DIF = EV_[0]-self.elpar[iel_][0]
        f_   = DIF*DIF
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
    def gSQUARE(self,nargout,*args):

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

