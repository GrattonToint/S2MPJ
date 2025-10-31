from s2mpjlib import *
class  LOTSCHD(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    A simple quadratic program inspired by the economic lot scheduling
#    problem.
# 
#    Source:
#    an exercize for L. Watson course on LANCELOT in the Spring 1993.
# 
#    SIF input: T. Kuan, Virginia Tech., Spring 1993.
# 
#    classification = "C-CQLR2-AN-12-7"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'LOTSCHD'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['X1'] = 1.502
        v_['X2'] = 1.126
        v_['X3'] = 0.815
        v_['X4'] = 1.268
        v_['X5'] = 1.502
        v_['X6'] = 0.740
        v_['A1'] = 1.8
        v_['A2'] = 3.2
        v_['A3'] = 6.1
        v_['A4'] = 3.2
        v_['A5'] = 1.8
        v_['A6'] = 7.4
        v_['C1'] = 11.0
        v_['C2'] = 3.0
        v_['C3'] = 20.0
        v_['C4'] = 17.0
        v_['C5'] = 9.0
        v_['C6'] = 20.0
        v_['C7'] = 126.1
        v_['1'] = 1
        v_['2'] = 2
        v_['4'] = 4
        v_['5'] = 5
        v_['6'] = 6
        v_['7'] = 7
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['6'])+1):
            [iv,ix_,_] = s2mpj_ii('T'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'T'+str(I))
            [iv,ix_,_] = s2mpj_ii('U'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'U'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['1']),int(v_['6'])+1):
            [ig,ig_,_] = s2mpj_ii('OBJ'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['T'+str(I)]])
            valA = np.append(valA,float(v_['X'+str(I)]))
            [ig,ig_,_] = s2mpj_ii('CONS7',ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'CONS7')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['T'+str(I)]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['U'+str(I)]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('CONS'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'CONS'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['T'+str(I)]])
            valA = np.append(valA,float(v_['A'+str(I)]))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['U'+str(I)]])
            valA = np.append(valA,float(-1.0))
        for I in range(int(v_['2']),int(v_['4'])+1):
            [ig,ig_,_] = s2mpj_ii('CONS'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'CONS'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['T'+str(I)]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['U'+str(I)]])
            valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('CONS2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'CONS2')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['T3']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['U3']])
        valA = np.append(valA,float(-1.0))
        for I in range(int(v_['1']),int(v_['2'])+1):
            [ig,ig_,_] = s2mpj_ii('CONS3',ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'CONS3')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['T'+str(I)]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['U'+str(I)]])
            valA = np.append(valA,float(-1.0))
        for I in range(int(v_['4']),int(v_['6'])+1):
            [ig,ig_,_] = s2mpj_ii('CONS3',ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'CONS3')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['T'+str(I)]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['U'+str(I)]])
            valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('CONS4',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'CONS4')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['T1']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['U1']])
        valA = np.append(valA,float(-1.0))
        for I in range(int(v_['5']),int(v_['6'])+1):
            [ig,ig_,_] = s2mpj_ii('CONS4',ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'CONS4')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['T'+str(I)]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['U'+str(I)]])
            valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('CONS5',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'CONS5')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['T1']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['U1']])
        valA = np.append(valA,float(-1.0))
        for I in range(int(v_['1']),int(v_['5'])+1):
            [ig,ig_,_] = s2mpj_ii('CONS6',ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'CONS6')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['T'+str(I)]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['U'+str(I)]])
            valA = np.append(valA,float(-1.0))
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
        for I in range(int(v_['1']),int(v_['7'])+1):
            self.gconst = arrset(self.gconst,ig_['CONS'+str(I)],float(v_['C'+str(I)]))
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
        for I in range(int(v_['1']),int(v_['6'])+1):
            ig = ig_['OBJ'+str(I)]
            self.grftype = arrset(self.grftype,ig,'gSQUARE')
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
        self.lincons   = np.arange(len(self.congrps))
        self.pbclass   = "C-CQLR2-AN-12-7"
        self.x0        = np.zeros((self.n,1))
        self.objderlvl = 2
        self.conderlvl = [2]


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

