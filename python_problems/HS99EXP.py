from s2mpjlib import *
class  HS99EXP(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS99EXP
#    *********
# 
#    Source: an expanded form of problem 99 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: Ph. Toint, April 1991.
# 
#    classification = "C-COOR2-AN-31-21"
# 
#    Problem data
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS99EXP'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['T1'] = 0.0
        v_['T2'] = 25.0
        v_['T3'] = 50.0
        v_['T4'] = 100.0
        v_['T5'] = 150.0
        v_['T6'] = 200.0
        v_['T7'] = 290.0
        v_['T8'] = 380.0
        v_['A1'] = 0.0
        v_['A2'] = 50.0
        v_['A3'] = 50.0
        v_['A4'] = 75.0
        v_['A5'] = 75.0
        v_['A6'] = 75.0
        v_['A7'] = 100.0
        v_['A8'] = 100.0
        v_['B'] = 32.0
        v_['1'] = 1
        v_['2'] = 2
        v_['7'] = 7
        v_['8'] = 8
        for I in range(int(v_['2']),int(v_['8'])+1):
            v_['I-1'] = -1+I
            v_['DT'+str(I)] = v_['T'+str(I)]-v_['T'+str(int(v_['I-1']))]
            v_['DTISQ'] = v_['DT'+str(I)]*v_['DT'+str(I)]
            v_['DT'+str(I)] = 0.5*v_['DTISQ']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['7'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I))
            [iv,ix_,_] = s2mpj_ii('R'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'R'+str(I))
            [iv,ix_,_] = s2mpj_ii('Q'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'Q'+str(I))
            [iv,ix_,_] = s2mpj_ii('S'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'S'+str(I))
        [iv,ix_,_] = s2mpj_ii('R'+str(int(v_['8'])),ix_)
        self.xnames=arrset(self.xnames,iv,'R'+str(int(v_['8'])))
        [iv,ix_,_] = s2mpj_ii('Q'+str(int(v_['8'])),ix_)
        self.xnames=arrset(self.xnames,iv,'Q'+str(int(v_['8'])))
        [iv,ix_,_] = s2mpj_ii('S'+str(int(v_['8'])),ix_)
        self.xnames=arrset(self.xnames,iv,'S'+str(int(v_['8'])))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['2']),int(v_['8'])+1):
            v_['I-1'] = -1+I
            [ig,ig_,_] = s2mpj_ii('R'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'R'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['R'+str(I)]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['R'+str(int(v_['I-1']))]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('Q'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'Q'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['Q'+str(I)]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['Q'+str(int(v_['I-1']))]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['S'+str(int(v_['I-1']))]])
            valA = np.append(valA,float(v_['DT'+str(I)]))
            [ig,ig_,_] = s2mpj_ii('S'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'S'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['S'+str(I)]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['S'+str(int(v_['I-1']))]])
            valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['R'+str(int(v_['8']))]])
        valA = np.append(valA,float(1.0))
        self.gscale = arrset(self.gscale,ig,float(-1.0))
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
        for I in range(int(v_['2']),int(v_['7'])+1):
            v_['RHS'] = v_['DT'+str(I)]*v_['B']
            self.gconst = arrset(self.gconst,ig_['Q'+str(I)],float(v_['RHS']))
            v_['RHS'] = v_['DT'+str(I)]*v_['B']
            self.gconst = arrset(self.gconst,ig_['S'+str(I)],float(v_['RHS']))
        self.gconst = arrset(self.gconst,ig_['Q'+str(int(v_['8']))],float(100000.0))
        self.gconst = arrset(self.gconst,ig_['S'+str(int(v_['8']))],float(1000.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        self.xlower[ix_['R'+str(int(v_['1']))]] = 0.0
        self.xupper[ix_['R'+str(int(v_['1']))]] = 0.0
        self.xlower[ix_['Q'+str(int(v_['1']))]] = 0.0
        self.xupper[ix_['Q'+str(int(v_['1']))]] = 0.0
        self.xlower[ix_['S'+str(int(v_['1']))]] = 0.0
        self.xupper[ix_['S'+str(int(v_['1']))]] = 0.0
        for I in range(int(v_['1']),int(v_['7'])+1):
            self.xlower[ix_['X'+str(I)]] = 0.0
            self.xupper[ix_['X'+str(I)]] = 1.58
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        for I in range(int(v_['1']),int(v_['7'])+1):
            self.x0[ix_['X'+str(I)]] = float(0.5)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSN', iet_)
        elftv = loaset(elftv,it,0,'X')
        [it,iet_,_] = s2mpj_ii( 'eCS', iet_)
        elftv = loaset(elftv,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['1']),int(v_['7'])+1):
            ename = 'SNX'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSN')
            ielftype = arrset(ielftype,ie,iet_["eSN"])
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'CSX'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eCS')
            ielftype = arrset(ielftype,ie,iet_["eCS"])
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
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
        ig = ig_['OBJ']
        self.grftype = arrset(self.grftype,ig,'gL2')
        for I in range(int(v_['2']),int(v_['8'])+1):
            v_['I-1'] = -1+I
            v_['W'] = v_['A'+str(I)]*v_['DT'+str(I)]
            ig = ig_['R'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['CSX'+str(int(v_['I-1']))])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['W']))
            ig = ig_['S'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['SNX'+str(int(v_['I-1']))])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['W']))
            v_['W'] = v_['A'+str(I)]*v_['DT'+str(I)]
            ig = ig_['Q'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['SNX'+str(int(v_['I-1']))])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['W']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               -831079892.0
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
        self.pbclass   = "C-COOR2-AN-31-21"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSN(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        SNX = np.sin(EV_[0,0])
        f_   = SNX
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = np.cos(EV_[0,0])
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = -SNX
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eCS(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        CSX = np.cos(EV_[0,0])
        f_   = CSX
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -np.sin(EV_[0,0])
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = -CSX
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

