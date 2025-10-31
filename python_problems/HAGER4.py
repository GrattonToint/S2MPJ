from s2mpjlib import *
class  HAGER4(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HAGER4
#    *********
# 
#    A nonlinear optimal control problem, by W. Hager.
# 
#    NOTE: The solution for x given in the article below by Hager has
#    a typo. On the interval [1/2, 1], x(t) = (exp(2t) + exp(t))/d. In
#    other words, the minus sign in the article should be a plus sign.
# 
#    Source: problem P4 in
#    W.W. Hager,
#    "Multiplier Methods for Nonlinear Optimal Control",
#    SIAM J. on Numercal Analysis 27(4): 1061-1080, 1990.
# 
#    SIF input: Ph. Toint, April 1991.
# 
#    classification = "C-COLR2-AN-V-V"
# 
#    Number of discretized points in [0,1]
# 
#           Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER     original value
# IE N                   50             $-PARAMETER
# IE N                   100            $-PARAMETER
# IE N                   500            $-PARAMETER
# IE N                   1000           $-PARAMETER
# IE N                   2500           $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HAGER4'

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
# IE N                   5000           $-PARAMETER
        v_['1/H'] = float(v_['N'])
        v_['H'] = 1.0/v_['1/H']
        v_['H/2'] = 0.5*v_['H']
        v_['1/H-1'] = -1.0+v_['1/H']
        v_['-1/H'] = -1.0*v_['1/H']
        v_['1/HSQ'] = v_['1/H']*v_['1/H']
        v_['1/2HSQ'] = 0.5*v_['1/HSQ']
        v_['0'] = 0
        v_['1'] = 1
        for I in range(int(v_['0']),int(v_['N'])+1):
            v_['RI'] = float(I)
            v_['T'+str(I)] = v_['RI']*v_['H']
            v_['-2TI'] = -2.0*v_['T'+str(I)]
            v_['Z'+str(I)] = np.exp(v_['-2TI'])
        for I in range(int(v_['0']),int(v_['1'])+1):
            v_['A'+str(I)] = -0.5*v_['Z'+str(I)]
            v_['TI+1/2'] = 0.5+v_['T'+str(I)]
            v_['B'+str(I)] = v_['A'+str(I)]*v_['TI+1/2']
            v_['TISQ'] = v_['T'+str(I)]*v_['T'+str(I)]
            v_['TIETC'] = v_['TISQ']+v_['TI+1/2']
            v_['C'+str(I)] = v_['A'+str(I)]*v_['TIETC']
        v_['DA'] = v_['A'+str(int(v_['1']))]-v_['A'+str(int(v_['0']))]
        v_['SCDA'] = 0.5*v_['DA']
        v_['DB'] = v_['B'+str(int(v_['1']))]-v_['B'+str(int(v_['0']))]
        v_['SCDB'] = v_['DB']*v_['1/H']
        v_['DC'] = v_['C'+str(int(v_['1']))]-v_['C'+str(int(v_['0']))]
        v_['SCDC'] = v_['DC']*v_['1/2HSQ']
        v_['E'] = np.exp(1.0)
        v_['3E'] = 3.0*v_['E']
        v_['1+3E'] = 1.0+v_['3E']
        v_['1-E'] = 1.0-v_['E']
        v_['2-2E'] = 2.0*v_['1-E']
        v_['XX0'] = v_['1+3E']/v_['2-2E']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['0']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I))
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('U'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'U'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['I-1'] = -1+I
            [ig,ig_,_] = s2mpj_ii('S'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'S'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(I)]])
            valA = np.append(valA,float(v_['1/H-1']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['I-1']))]])
            valA = np.append(valA,float(v_['-1/H']))
            v_['ETI'] = np.exp(v_['T'+str(I)])
            v_['-ETI'] = -1.0*v_['ETI']
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['U'+str(I)]])
            valA = np.append(valA,float(v_['-ETI']))
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
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        self.xlower[ix_['X'+str(int(v_['0']))]] = v_['XX0']
        self.xupper[ix_['X'+str(int(v_['0']))]] = v_['XX0']
        for I in range(int(v_['1']),int(v_['N'])+1):
            self.xupper[ix_['U'+str(I)]] = 1.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        self.x0[ix_['X'+str(int(v_['0']))]] = float(v_['XX0'])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eELT', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        elftp = []
        elftp = loaset(elftp,it,0,'D')
        elftp = loaset(elftp,it,1,'E')
        elftp = loaset(elftp,it,2,'F')
        [it,iet_,_] = s2mpj_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['I-1'] = -1+I
            ename = 'EL'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eELT')
            ielftype = arrset(ielftype,ie,iet_["eELT"])
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I-1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='Y')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            v_['DD'] = v_['SCDA']*v_['Z'+str(int(v_['I-1']))]
            v_['EE'] = v_['SCDB']*v_['Z'+str(int(v_['I-1']))]
            v_['FF'] = v_['SCDC']*v_['Z'+str(int(v_['I-1']))]
            posep = np.where(elftp[ielftype[ie]]=='D')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['DD']))
            posep = np.where(elftp[ielftype[ie]]=='E')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['EE']))
            posep = np.where(elftp[ielftype[ie]]=='F')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['FF']))
            ename = 'U'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSQ')
            ielftype = arrset(ielftype,ie,iet_["eSQ"])
            vname = 'U'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            ig = ig_['OBJ']
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['EL'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['U'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['H/2']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN(10)           2.833914199
# LO SOLTN(50)           2.799810928
# LO SOLTN(100)          2.796761851
# LO SOLTN(500)          2.794513229
# LO SOLTN(1000)         2.794244187
# LO SOLTN(5000)         ???
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-COLR2-AN-V-V"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eELT(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_    = (
              self.elpar[iel_][0]*EV_[0]*EV_[0]+self.elpar[iel_][1]*EV_[0]*(EV_[1]-EV_[0])+self.elpar[iel_][2]*(EV_[1]-EV_[0])**2)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = (2.0*self.elpar[iel_][0]*EV_[0]+self.elpar[iel_][1]*(EV_[1]-2.0*EV_[0])-
                 2.0*self.elpar[iel_][2]*(EV_[1]-EV_[0]))
            g_[1] = self.elpar[iel_][1]*EV_[0]+2.0*self.elpar[iel_][2]*(EV_[1]-EV_[0])
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 2.0*(self.elpar[iel_][0]-self.elpar[iel_][1]+self.elpar[iel_][2])
                H_[0,1] = self.elpar[iel_][1]-2.0*self.elpar[iel_][2]
                H_[1,0] = H_[0,1]
                H_[1,1] = 2.0*self.elpar[iel_][2]
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
        f_   = EV_[0]*EV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[0]+EV_[0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

