from s2mpjlib import *
class  TORSION6(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : TORSION6
#    *********
# 
#    The quadratic elastic torsion problem
# 
#    The problem comes from the obstacle problem on a square.
# 
#    The square is discretized into (px-1)(py-1) little squares. The
#    heights of the considered surface above the corners of these little
#    squares are the problem variables,  There are px**2 of them.
# 
#    The dimension of the problem is specified by Q, which is half the
#    number discretization points along one of the coordinate
#    direction.
#    Since the number of variables is P**2, it is given by 4Q**2
# 
#    Source: problem (c=20, starting point Z = origin) in
#    J. More' and G. Toraldo,
#    "On the Solution of Large Quadratic-Programming Problems with Bound
#    Constraints", 
#    SIAM J. on Optimization, vol 1(1), pp. 93-113, 1991.
# 
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "C-CQBR2-MY-V-0"
# 
#    Q is half the number of discretized points along the X axis
# 
#           Alternative values for the SIF file parameters:
# IE Q                   2              $-PARAMETER n= 16
# IE Q                   5              $-PARAMETER n= 100     original value
# IE Q                   11             $-PARAMETER n= 484
# IE Q                   16             $-PARAMETER n= 1024
# IE Q                   37             $-PARAMETER n= 5476
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'TORSION6'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['Q'] = int(2);  #  SIF file default value
        else:
            v_['Q'] = int(args[0])
# IE Q                   50             $-PARAMETER n= 10000
# IE Q                   61             $-PARAMETER n= 14884
        if nargin<2:
            v_['C'] = float(20.0);  #  SIF file default value
        else:
            v_['C'] = float(args[1])
        v_['Q+1'] = 1+v_['Q']
        v_['P'] = v_['Q']+v_['Q']
        v_['P-1'] = -1+v_['P']
        v_['1/H'] = float(v_['P-1'])
        v_['H'] = 1.0/v_['1/H']
        v_['H2'] = v_['H']*v_['H']
        v_['C0'] = v_['H2']*v_['C']
        v_['LC'] = -1.0*v_['C0']
        v_['1'] = 1
        v_['2'] = 2
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
        for I in range(int(v_['2']),int(v_['P-1'])+1):
            for J in range(int(v_['2']),int(v_['P-1'])+1):
                [ig,ig_,_] = s2mpj_ii('G'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(I)+','+str(J)]])
                valA = np.append(valA,float(v_['LC']))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        selfnob      = ngrp
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        for J in range(int(v_['1']),int(v_['P'])+1):
            self.xlower[ix_['X'+str(int(v_['1']))+','+str(J)]] = 0.0
            self.xupper[ix_['X'+str(int(v_['1']))+','+str(J)]] = 0.0
            self.xlower[ix_['X'+str(int(v_['P']))+','+str(J)]] = 0.0
            self.xupper[ix_['X'+str(int(v_['P']))+','+str(J)]] = 0.0
        for I in range(int(v_['2']),int(v_['P-1'])+1):
            self.xlower[ix_['X'+str(I)+','+str(int(v_['P']))]] = 0.0
            self.xupper[ix_['X'+str(I)+','+str(int(v_['P']))]] = 0.0
            self.xlower[ix_['X'+str(I)+','+str(int(v_['1']))]] = 0.0
            self.xupper[ix_['X'+str(I)+','+str(int(v_['1']))]] = 0.0
        for I in range(int(v_['2']),int(v_['Q'])+1):
            for J in range(int(v_['2']),int(I)+1):
                v_['J-1'] = -1+J
                v_['RJ-1'] = float(v_['J-1'])
                v_['UPPL'] = v_['RJ-1']*v_['H']
                v_['LOWL'] = -1.0*v_['UPPL']
                self.xlower[ix_['X'+str(I)+','+str(J)]] = v_['LOWL']
                self.xupper[ix_['X'+str(I)+','+str(J)]] = v_['UPPL']
            v_['MI'] = -1*I
            v_['P-I'] = v_['P']+v_['MI']
            v_['I-1'] = -1+I
            v_['RI-1'] = float(v_['I-1'])
            v_['UPPM'] = v_['RI-1']*v_['H']
            v_['LOWM'] = -1.0*v_['UPPM']
            v_['P-I+1'] = 1+v_['P-I']
            for J in range(int(I),int(v_['P-I+1'])+1):
                self.xlower[ix_['X'+str(I)+','+str(J)]] = v_['LOWM']
                self.xupper[ix_['X'+str(I)+','+str(J)]] = v_['UPPM']
            for J in range(int(v_['P-I+1']),int(v_['P-1'])+1):
                v_['MJ'] = -1*J
                v_['P-J'] = v_['P']+v_['MJ']
                v_['RP-J'] = float(v_['P-J'])
                v_['UPPR'] = v_['RP-J']*v_['H']
                v_['LOWR'] = -1.0*v_['UPPR']
                self.xlower[ix_['X'+str(I)+','+str(J)]] = v_['LOWR']
                self.xupper[ix_['X'+str(I)+','+str(J)]] = v_['UPPR']
        for I in range(int(v_['Q+1']),int(v_['P-1'])+1):
            v_['MI'] = -1*I
            v_['P-I'] = v_['P']+v_['MI']
            v_['P-I+1'] = 1+v_['P-I']
            for J in range(int(v_['2']),int(v_['P-I+1'])+1):
                v_['J-1'] = -1+J
                v_['RJ-1'] = float(v_['J-1'])
                v_['UPPL'] = v_['RJ-1']*v_['H']
                v_['LOWL'] = -1.0*v_['UPPL']
                self.xlower[ix_['X'+str(I)+','+str(J)]] = v_['LOWL']
                self.xupper[ix_['X'+str(I)+','+str(J)]] = v_['UPPL']
            v_['RP-I'] = float(v_['P-I'])
            v_['UPPM'] = v_['RP-I']*v_['H']
            v_['LOWM'] = -1.0*v_['UPPM']
            for J in range(int(v_['P-I+1']),int(I)+1):
                self.xlower[ix_['X'+str(I)+','+str(J)]] = v_['LOWM']
                self.xupper[ix_['X'+str(I)+','+str(J)]] = v_['UPPM']
            for J in range(int(I),int(v_['P-1'])+1):
                v_['MJ'] = -1*J
                v_['P-J'] = v_['P']+v_['MJ']
                v_['RP-J'] = float(v_['P-J'])
                v_['UPPR'] = v_['RP-J']*v_['H']
                v_['LOWR'] = -1.0*v_['UPPR']
                self.xlower[ix_['X'+str(I)+','+str(J)]] = v_['LOWR']
                self.xupper[ix_['X'+str(I)+','+str(J)]] = v_['UPPR']
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
        for I in range(int(v_['2']),int(v_['P-1'])+1):
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            for J in range(int(v_['2']),int(v_['P-1'])+1):
                v_['J-1'] = -1+J
                v_['J+1'] = 1+J
                ename = 'A'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eISQ')
                ielftype = arrset(ielftype,ie,iet_["eISQ"])
                self.x0 = np.zeros((self.n,1))
                vname = 'X'+str(int(v_['I+1']))+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                ename = 'B'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eISQ')
                ielftype = arrset(ielftype,ie,iet_["eISQ"])
                vname = 'X'+str(I)+','+str(int(v_['J+1']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                ename = 'C'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eISQ')
                ielftype = arrset(ielftype,ie,iet_["eISQ"])
                vname = 'X'+str(int(v_['I-1']))+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                ename = 'D'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eISQ')
                ielftype = arrset(ielftype,ie,iet_["eISQ"])
                vname = 'X'+str(I)+','+str(int(v_['J-1']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['2']),int(v_['P-1'])+1):
            for J in range(int(v_['2']),int(v_['P-1'])+1):
                ig = ig_['G'+str(I)+','+str(J)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['A'+str(I)+','+str(J)])
                self.grelw = loaset(self.grelw,ig,posel,float(0.25))
                posel = posel+1
                self.grelt = loaset(self.grelt,ig,posel,ie_['B'+str(I)+','+str(J)])
                self.grelw = loaset(self.grelw,ig,posel,float(0.25))
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['C'+str(I)+','+str(J)])
                self.grelw = loaset(self.grelw,ig,posel,float(0.25))
                posel = posel+1
                self.grelt = loaset(self.grelt,ig,posel,ie_['D'+str(I)+','+str(J)])
                self.grelw = loaset(self.grelw,ig,posel,float(0.25))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN(2)            -2.7407407353
# LO SOLTN(5)            -2.8971193358
# LO SOLTN(11)           -2.8847068155
# LO SOLTN(16)           ???
# LO SOLTN(37)           ???
# LO SOLTN(50)           ???
# LO SOLTN(61)           ???
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass   = "C-CQBR2-MY-V-0"
        self.x0        = np.zeros((self.n,1))
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

