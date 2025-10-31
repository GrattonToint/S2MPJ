from s2mpjlib import *
class  CLPLATEB(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CLPLATEB
#    *********
# 
#    The clamped plate problem (Strang, Nocedal, Dax)
#    The problem comes from the discretization the following problem
#    in mechanics:  a plate is clamped on one edge and loaded on the
#    opposite side.  The plate is the unit square.
# 
#    The plate is clamped on its lower edge, by fixing the
#    corresponding variables to zero.
# 
#    In this version of the problem, the weight WGHT is distributed
#    equally along the upper edge, introducing a symmetry with respect
#    to the vertical axis.
# 
#    Source:
#    J. Nocedal,
#    "Solving large nonlinear systems of equations arising in mechanics",
#    Proceedings of the Cocoyoc Numerical Analysis Conference, Mexico,
#    pp. 132-141, 1981.
# 
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "C-COXR2-MN-V-0"
# 
#    P is the number of points in one side of the unit square
#    The number of variables is P*P, of which (P-1)*(P-1) are free.
# 
#           Alternative values for the SIF file parameters:
# IE P                   4              $-PARAMETER n = 16
# IE P                   7              $-PARAMETER n = 49    original value
# IE P                   10             $-PARAMETER n = 100
# IE P                   23             $-PARAMETER n = 529
# IE P                   32             $-PARAMETER n = 1024
# IE P                   71             $-PARAMETER n = 5041
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'CLPLATEB'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['P'] = int(4);  #  SIF file default value
        else:
            v_['P'] = int(args[0])
        v_['WGHT'] = -0.1
        v_['1'] = 1
        v_['2'] = 2
        v_['P-1'] = -1+v_['P']
        v_['RP-1'] = float(v_['P-1'])
        v_['1/P-1'] = 1.0/v_['RP-1']
        v_['DISW'] = v_['WGHT']*v_['1/P-1']
        v_['P2'] = v_['P']*v_['P']
        v_['RP2'] = float(v_['P2'])
        v_['HP2'] = 0.5*v_['RP2']
        v_['1/HP2'] = 1.0/v_['HP2']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for J in range(int(v_['1']),int(v_['P'])+1):
            for I in range(int(v_['1']),int(v_['P'])+1):
                [iv,ix_,_] = s2mpj_ii('X'+str(I)+','+str(J),ix_)
                self.xnames=arrset(self.xnames,iv,'X'+str(I)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['2']),int(v_['P'])+1):
            v_['I-1'] = -1+I
            for J in range(int(v_['2']),int(v_['P'])+1):
                v_['J-1'] = -1+J
                [ig,ig_,_] = s2mpj_ii('A'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
                self.gscale = arrset(self.gscale,ig,float(2.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(I)+','+str(J)]])
                valA = np.append(valA,float(1.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(I)+','+str(int(v_['J-1']))]])
                valA = np.append(valA,float(-1.0))
                [ig,ig_,_] = s2mpj_ii('B'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
                self.gscale = arrset(self.gscale,ig,float(2.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(I)+','+str(J)]])
                valA = np.append(valA,float(1.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(int(v_['I-1']))+','+str(J)]])
                valA = np.append(valA,float(-1.0))
                [ig,ig_,_] = s2mpj_ii('C'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
                self.gscale = arrset(self.gscale,ig,float(v_['1/HP2']))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(I)+','+str(J)]])
                valA = np.append(valA,float(1.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(I)+','+str(int(v_['J-1']))]])
                valA = np.append(valA,float(-1.0))
                [ig,ig_,_] = s2mpj_ii('D'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
                self.gscale = arrset(self.gscale,ig,float(v_['1/HP2']))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(I)+','+str(J)]])
                valA = np.append(valA,float(1.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(int(v_['I-1']))+','+str(J)]])
                valA = np.append(valA,float(-1.0))
        for J in range(int(v_['1']),int(v_['P'])+1):
            [ig,ig_,_] = s2mpj_ii('W',ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['P']))+','+str(J)]])
            valA = np.append(valA,float(v_['DISW']))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        for J in range(int(v_['1']),int(v_['P'])+1):
            self.xlower[ix_['X'+str(int(v_['1']))+','+str(J)]] = 0.0
            self.xupper[ix_['X'+str(int(v_['1']))+','+str(J)]] = 0.0
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(0.0))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gL2',igt_)
        [it,igt_,_] = s2mpj_ii('gL4',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['2']),int(v_['P'])+1):
            for J in range(int(v_['2']),int(v_['P'])+1):
                ig = ig_['A'+str(I)+','+str(J)]
                self.grftype = arrset(self.grftype,ig,'gL2')
                ig = ig_['B'+str(I)+','+str(J)]
                self.grftype = arrset(self.grftype,ig,'gL2')
                ig = ig_['C'+str(I)+','+str(J)]
                self.grftype = arrset(self.grftype,ig,'gL4')
                ig = ig_['D'+str(I)+','+str(J)]
                self.grftype = arrset(self.grftype,ig,'gL4')
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 0.0
#    Solution
# LO SOLTN(4)            -9.3705D-03
# LO SOLTN(7)            -6.9193D-03
# LO SOLTN(10)           -6.2008D-03
# LO SOLTN(23)           -5.4274D-03
# LO SOLTN(32)           -5.2835D-03
# LO SOLTN(71)           -5.0948D-03
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass   = "C-COXR2-MN-V-0"
        self.objderlvl = 2

# ********************
#  SET UP THE GROUPS *
#  ROUTINE           *
# ********************

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

    @staticmethod
    def gL4(self,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= GVAR_**4
        if nargout>1:
            g_ = 4.0*GVAR_**3
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 12.0*GVAR_**2
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

