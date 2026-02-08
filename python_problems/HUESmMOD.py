from s2mpjlib import *
class  HUESmMOD(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HUESmMOD
#    *********
# 
#    Source: An inverse problem from astronomy,
#    reformulated as a convex quadratic program by
#    S. P. Hestis, SIAM Review 34 (1992) pp. 642-647.
# 
#    SIF input: Nick Gould, January 1993.
#    improvements by: Ruediger Franke (Ruediger.Franke@RZ.TU-Ilmenau.DE)
# 
#    classification = "C-CQLR2-MN-V-V"
# 
#    Number of variables
# 
#           Alternative values for the SIF file parameters:
# IE K                   10             $-PARAMETER
# IE K                   100            $-PARAMETER
# IE K                   1000           $-PARAMETER    original value
# IE K                   5000           $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HUESmMOD'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['K'] = int(10);  #  SIF file default value
        else:
            v_['K'] = int(args[0])
# IE K                   10000          $-PARAMETER
        v_['1'] = 1
        v_['RANGE'] = 1.0
        v_['3.0'] = 3.0
        v_['5.0'] = 5.0
        v_['RK'] = float(v_['K'])
        v_['DELTAX'] = v_['RANGE']/v_['RK']
        v_['DELTAX2'] = v_['DELTAX']*v_['DELTAX']
        v_['DELTAX3'] = v_['DELTAX2']*v_['DELTAX']
        v_['DELTAX5'] = v_['DELTAX3']*v_['DELTAX2']
        v_['DELTAX3/3'] = v_['DELTAX3']/v_['3.0']
        v_['DELTAX5/5'] = v_['DELTAX5']/v_['5.0']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['K'])+1):
            [iv,ix_,_] = s2mpj_ii('M'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'M'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['1']),int(v_['K'])+1):
            [ig,ig_,_] = s2mpj_ii('OBJ'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['1']),int(v_['K'])+1):
            v_['I-1'] = -1+I
            v_['RI'] = float(I)
            v_['RI2'] = v_['RI']*v_['RI']
            v_['RI3'] = v_['RI2']*v_['RI']
            v_['RI5'] = v_['RI3']*v_['RI2']
            v_['RI-1'] = float(v_['I-1'])
            v_['RI-12'] = v_['RI-1']*v_['RI-1']
            v_['RI-13'] = v_['RI-12']*v_['RI-1']
            v_['RI-15'] = v_['RI-13']*v_['RI-12']
            v_['DIFF3'] = v_['RI3']-v_['RI-13']
            v_['DIFF5'] = v_['RI5']-v_['RI-15']
            v_['COEFF1'] = v_['DIFF3']*v_['DELTAX3/3']
            v_['COEFF2'] = v_['DIFF5']*v_['DELTAX5/5']
            [ig,ig_,_] = s2mpj_ii('E1',ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'E1')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['M'+str(I)]])
            valA = np.append(valA,float(v_['COEFF1']))
            [ig,ig_,_] = s2mpj_ii('E2',ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'E2')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['M'+str(I)]])
            valA = np.append(valA,float(v_['COEFF2']))
        v_['RK'] = float(v_['K'])
        v_['1/RK'] = 1.0/v_['RK']
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
        self.gconst = arrset(self.gconst,ig_['E1'],float(1835.2))
        self.gconst = arrset(self.gconst,ig_['E2'],float(909.8))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(1.0))
        self.y0 = np.full((self.m,1),float(1.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'U1')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['1']),int(v_['K'])+1):
            ename = 'E'+str(I)
            [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
            if newelt:
                self.elftype = arrset(self.elftype,ie,'eSQ')
                ielftype = arrset(ielftype,ie,iet_['eSQ'])
            vname = 'M'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(1.0))
            posev = np.where(elftv[ielftype[ie]]=='U1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['K'])+1):
            ig = ig_['OBJ'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['1/RK']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 0.0
#    Solution
# LO SOLTN               0.0
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CQLR2-MN-V-V"
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

