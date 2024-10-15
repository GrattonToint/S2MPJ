from s2mpjlib import *
class  LISWET5(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : LISWET5
#    *********
# 
#    A k-convex approximation problem posed as a 
#    convex quadratic problem, with variable dimensions.
# 
#    Formulation:
#    -----------
# 
#                 n+k             2
#    minimize 1/2 sum ( x  - c  )
#                 i=1    i    i
# 
#    subject to
# 
#                  k              k-i
#                 sum ( k ) ( -1 )    x     > 0
#                 i=0 ( i )            j+i  = 
# 
#    where c  = g( t ) + small perturbation, t  = (i-1)/(n+k-1)
#           i       i                         i 
# 
#    Case 5: g(t) = exp(t)
# 
#    NB. Perturbations are not random as Li and Swetits's 
#        random number generator is undefined.
# 
#    Source:
#    W. Li and J. Swetits,
#    "A Newton method for convex regression, data smoothing and
#    quadratic programming with bounded constraints",
#    SIAM J. Optimization 3 (3) pp 466-488, 1993.
# 
#    SIF input: Nick Gould, August 1994.
# 
#    classification = "C-QLR2-AN-V-V"
# 
#           Alternative values for the SIF file parameters:
# IE N                   100            $-PARAMETER 103 variables original value 
# IE K                   3              $-PARAMETER original value
# 
# IE N                   100            $-PARAMETER 104 variables    
# IE K                   4              $-PARAMETER
# 
# IE N                   100            $-PARAMETER 105 variables    
# IE K                   5              $-PARAMETER
# 
# IE N                   100            $-PARAMETER 106 variables    
# IE K                   6              $-PARAMETER
# 
# IE N                   400            $-PARAMETER 402 variables    
# IE K                   2              $-PARAMETER
# 
# IE N                   400            $-PARAMETER 403 variables    
# IE K                   3              $-PARAMETER
# 
# IE N                   2000           $-PARAMETER 2001 variables    
# IE K                   1              $-PARAMETER
# 
# IE N                   2000           $-PARAMETER 2002 variables    
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'LISWET5'

    def __init__(self, *args): 
        import numpy as np
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(50);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
# IE K                   2              $-PARAMETER
        if nargin<2:
            v_['K'] = int(3);  #  SIF file default value
        else:
            v_['K'] = int(args[1])
# IE N                   10000          $-PARAMETER 10001 variables    
# IE K                   1              $-PARAMETER
# IE N                   10000          $-PARAMETER 10002 variables    
# IE K                   2              $-PARAMETER
        v_['0'] = 0
        v_['1'] = 1
        v_['ONE'] = 1.0
        v_['HALF'] = 0.5
        v_['N+K'] = v_['N']+v_['K']
        v_['N+K-1'] = -1+v_['N+K']
        v_['RN+K-1'] = float(v_['N+K-1'])
        v_['CONST'] = 0.0
        v_['B'+str(int(v_['0']))] = v_['ONE']
        for I in range(int(v_['1']),int(v_['K'])+1):
            v_['I-1'] = -1+I
            v_['RI'] = float(I)
            v_['B'+str(I)] = v_['B'+str(int(v_['I-1']))]*v_['RI']
        v_['C'+str(int(v_['0']))] = v_['ONE']
        v_['PLUSMINUS'] = v_['ONE']
        for I in range(int(v_['1']),int(v_['K'])+1):
            v_['K-I'] = v_['K']-I
            v_['PLUSMINUS'] = -1.0*v_['PLUSMINUS']
            v_['C'+str(I)] = v_['B'+str(int(v_['K']))]/v_['B'+str(I)]
            v_['C'+str(I)] = v_['C'+str(I)]/v_['B'+str(int(v_['K-I']))]
            v_['C'+str(I)] = v_['C'+str(I)]*v_['PLUSMINUS']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N+K'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.A       = lil_matrix((1000000,1000000))
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames      = np.array([])
        self.cnames = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['N+K'])+1):
            v_['I-1'] = -1+I
            v_['RI'] = float(I)
            v_['RI-1'] = float(v_['I-1'])
            v_['TI'] = v_['RI-1']/v_['RN+K-1']
            v_['GT'] = np.exp(v_['TI'])
            v_['RANDOM'] = np.sin(v_['RI'])
            v_['RANDOM'] = 0.1*v_['RANDOM']
            v_['CI'] = v_['GT']+v_['RANDOM']
            v_['-CI'] = -1.0*v_['CI']
            v_['-CI*CI'] = v_['-CI']*v_['CI']
            v_['CONST'] = v_['CONST']+v_['-CI*CI']
            [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(I)]
            self.A[ig,iv] = float(v_['-CI'])+self.A[ig,iv]
        for J in range(int(v_['1']),int(v_['N'])+1):
            v_['J+K'] = J+v_['K']
            for I in range(int(v_['0']),int(v_['K'])+1):
                v_['J+K-I'] = v_['J+K']-I
                [ig,ig_,_] = s2mpj_ii('CON'+str(J),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'CON'+str(J))
                iv = ix_['X'+str(int(v_['J+K-I']))]
                self.A[ig,iv] = float(v_['C'+str(I)])+self.A[ig,iv]
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
        self.cnames= cnames[self.congrps]
        self.nob = ngrp-self.m
        self.objgrps = np.where(gtype=='<>')[0]
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        v_['CONST'] = v_['HALF']*v_['CONST']
        self.gconst = arrset(self.gconst,ig_['OBJ'],float(v_['CONST']))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        self.xlower = np.zeros((self.n,1))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['1']),int(v_['N+K'])+1):
            ename = 'XSQ'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSQ')
            ielftype = arrset(ielftype,ie,iet_["eSQ"])
            self.x0 = np.zeros((self.n,1))
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['N+K'])+1):
            ig = ig_['OBJ']
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['XSQ'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle+self.neq,self.m)] = np.zeros((self.nge,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        self.A.resize(ngrp,self.n)
        self.A     = self.A.tocsr()
        sA1,sA2    = self.A.shape
        self.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass = "C-QLR2-AN-V-V"
        self.x0        = np.zeros((self.n,1))
        self.objderlvl = 2
        self.conderlvl = [2]

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
        f_   = 5.0e-1*EV_[0]*EV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 1.0e+0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

