from s2mpjlib import *
class  LMINSURF(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : LMINSURF
#    *********
# 
#    The linear minimum surface problem.
# 
#    The problem comes from the discretization of the minimum surface
#    problem on the unit square: given a set of boundary conditions on
#    the four sides of the square, one must find the surface which
#    meets these boundary conditions and is of minimum area.
# 
#    The unit square is discretized into (p-1)**2 little squares. The
#    heights of the considered surface above the corners of these little
#    squares are the problem variables,  There are p**2 of them.
#    Given these heights, the area above a little square is
#    approximated by the
#      S(i,j) = sqrt( 1 + 0.5(p-1)**2 ( a(i,j)**2 + b(i,j)**2 ) ) / (p-1)**2
#    where
#      a(i,j) = x(i,j) - x(i+1,j+1)
#    and
#      b(i,j) = x(i+1,j) - x(i,j+1)
# 
#    In the Linear Mininum Surface, the boundary conditions are given
#    as the heights of a given plane above the square boundaries.  This
#    plane is specified by its height above the (0,0) point (H00 below),
#    and its slopes along the first and second coordinate
#    directions in the plane (these slopes are denoted SLOPEJ and SLOPEI below).
# 
#    Source:
#    A Griewank and Ph. Toint,
#    "Partitioned variable metric updates for large structured
#    optimization problems",
#    Numerische Mathematik 39:429-448, 1982.
# 
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "C-COXR2-MY-V-0"
# 
#    P is the number of points in one side of the unit square
# 
#           Alternative values for the SIF file parameters:
# IE P                   4              $-PARAMETER n = 16     original value
# IE P                   7              $-PARAMETER n = 49
# IE P                   8              $-PARAMETER n = 64
# IE P                   11             $-PARAMETER n = 121
# IE P                   31             $-PARAMETER n = 961
# IE P                   32             $-PARAMETER n = 1024
# IE P                   75             $-PARAMETER n = 5625
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'LMINSURF'

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
# IE P                   100            $-PARAMETER n = 10000
# IE P                   125            $-PARAMETER n = 15625
        v_['H00'] = 1.0
        v_['SLOPEJ'] = 4.0
        v_['SLOPEI'] = 8.0
        v_['TWOP'] = v_['P']+v_['P']
        v_['P-1'] = -1+v_['P']
        v_['PP-1'] = v_['P']*v_['P-1']
        v_['RP-1'] = float(v_['P-1'])
        v_['INVP-1'] = 1.0/v_['RP-1']
        v_['RP-1SQ'] = v_['INVP-1']*v_['INVP-1']
        v_['SCALE'] = 1.0/v_['RP-1SQ']
        v_['SQP-1'] = v_['RP-1']*v_['RP-1']
        v_['PARAM'] = 0.5*v_['SQP-1']
        v_['1'] = 1
        v_['2'] = 2
        v_['STON'] = v_['INVP-1']*v_['SLOPEI']
        v_['WTOE'] = v_['INVP-1']*v_['SLOPEJ']
        v_['H01'] = v_['H00']+v_['SLOPEJ']
        v_['H10'] = v_['H00']+v_['SLOPEI']
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
        for I in range(int(v_['1']),int(v_['P-1'])+1):
            for J in range(int(v_['1']),int(v_['P-1'])+1):
                [ig,ig_,_] = s2mpj_ii('S'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
                self.gscale = arrset(self.gscale,ig,float(v_['SCALE']))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        selfnob      = ngrp
        #%%%%%%%%%%%%%%%%%%  CONSTANTS %%%%%%%%%%%%%%%%%%%
        self.gconst = np.full((ngrp,1),-1.0)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        for J in range(int(v_['1']),int(v_['P'])+1):
            v_['J-1'] = -1+J
            v_['RJ-1'] = float(v_['J-1'])
            v_['TH'] = v_['RJ-1']*v_['WTOE']
            v_['TL'] = v_['TH']+v_['H00']
            v_['TU'] = v_['TH']+v_['H10']
            self.xlower[ix_['X'+str(int(v_['1']))+','+str(J)]] = v_['TL']
            self.xupper[ix_['X'+str(int(v_['1']))+','+str(J)]] = v_['TL']
            self.xlower[ix_['X'+str(int(v_['P']))+','+str(J)]] = v_['TU']
            self.xupper[ix_['X'+str(int(v_['P']))+','+str(J)]] = v_['TU']
        for I in range(int(v_['2']),int(v_['P-1'])+1):
            v_['I-1'] = -1+I
            v_['RI-1'] = float(v_['I-1'])
            v_['TV'] = v_['RI-1']*v_['STON']
            v_['TR'] = v_['TV']+v_['H00']
            v_['TL'] = v_['TV']+v_['H01']
            self.xlower[ix_['X'+str(I)+','+str(int(v_['P']))]] = v_['TL']
            self.xupper[ix_['X'+str(I)+','+str(int(v_['P']))]] = v_['TL']
            self.xlower[ix_['X'+str(I)+','+str(int(v_['1']))]] = v_['TR']
            self.xupper[ix_['X'+str(I)+','+str(int(v_['1']))]] = v_['TR']
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(0.0))
        for J in range(int(v_['1']),int(v_['P'])+1):
            v_['J-1'] = -1+J
            v_['RJ-1'] = float(v_['J-1'])
            v_['TH'] = v_['RJ-1']*v_['WTOE']
            v_['TL'] = v_['TH']+v_['H00']
            v_['TU'] = v_['TH']+v_['H10']
            self.x0[ix_['X'+str(int(v_['1']))+','+str(J)]] = float(v_['TL'])
            self.x0[ix_['X'+str(int(v_['P']))+','+str(J)]] = float(v_['TU'])
        for I in range(int(v_['2']),int(v_['P-1'])+1):
            v_['I-1'] = -1+I
            v_['RI-1'] = float(v_['I-1'])
            v_['TV'] = v_['RI-1']*v_['STON']
            v_['TR'] = v_['TV']+v_['H00']
            v_['TL'] = v_['TV']+v_['H01']
            self.x0[ix_['X'+str(I)+','+str(int(v_['P']))]] = float(v_['TL'])
            self.x0[ix_['X'+str(I)+','+str(int(v_['1']))]] = float(v_['TR'])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eISQ', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['1']),int(v_['P-1'])+1):
            v_['I+1'] = 1+I
            for J in range(int(v_['1']),int(v_['P-1'])+1):
                v_['J+1'] = 1+J
                ename = 'A'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eISQ')
                ielftype = arrset(ielftype,ie,iet_["eISQ"])
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
                posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'X'+str(int(v_['I+1']))+','+str(int(v_['J+1']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
                posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                ename = 'B'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eISQ')
                ielftype = arrset(ielftype,ie,iet_["eISQ"])
                vname = 'X'+str(int(v_['I+1']))+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
                posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'X'+str(I)+','+str(int(v_['J+1']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
                posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gSQROOT',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['P-1'])+1):
            for J in range(int(v_['1']),int(v_['P-1'])+1):
                ig = ig_['S'+str(I)+','+str(J)]
                self.grftype = arrset(self.grftype,ig,'gSQROOT')
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['A'+str(I)+','+str(J)])
                self.grelw = loaset(self.grelw,ig,posel,float(v_['PARAM']))
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['B'+str(I)+','+str(J)])
                self.grelw = loaset(self.grelw,ig,posel,float(v_['PARAM']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 0.0
#    Solution
# LO SOLTN               9.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass   = "C-COXR2-MY-V-0"
        self.objderlvl = 2

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eISQ(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((1,2))
        IV_ = np.zeros(1)
        U_[0,0] = U_[0,0]+1
        U_[0,1] = U_[0,1]-1
        IV_[0] = to_scalar(U_[0:1,:].dot(EV_))
        f_   = IV_[0]*IV_[0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = IV_[0]+IV_[0]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gSQROOT(self,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        SQRAL = np.sqrt(GVAR_)
        f_= SQRAL
        if nargout>1:
            g_ = 0.5e0/SQRAL
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = -0.25e0/(SQRAL*GVAR_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

