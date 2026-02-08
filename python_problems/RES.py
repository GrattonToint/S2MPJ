from s2mpjlib import *
class  RES(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : RES
#    *********
# 
#    Dassault France ressort (spring) problem
# 
#    SIF input:  A. R. Conn, June 1993.
# 
#    classification = "C-CNLR2-MN-20-14"
# 
# 
# 
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'RES'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [iv,ix_,_] = s2mpj_ii('L0',ix_)
        self.xnames=arrset(self.xnames,iv,'L0')
        [iv,ix_,_] = s2mpj_ii('N',ix_)
        self.xnames=arrset(self.xnames,iv,'N')
        [iv,ix_,_] = s2mpj_ii('F',ix_)
        self.xnames=arrset(self.xnames,iv,'F')
        [iv,ix_,_] = s2mpj_ii('K',ix_)
        self.xnames=arrset(self.xnames,iv,'K')
        [iv,ix_,_] = s2mpj_ii('LB',ix_)
        self.xnames=arrset(self.xnames,iv,'LB')
        [iv,ix_,_] = s2mpj_ii('L',ix_)
        self.xnames=arrset(self.xnames,iv,'L')
        [iv,ix_,_] = s2mpj_ii('DE',ix_)
        self.xnames=arrset(self.xnames,iv,'DE')
        [iv,ix_,_] = s2mpj_ii('DI',ix_)
        self.xnames=arrset(self.xnames,iv,'DI')
        [iv,ix_,_] = s2mpj_ii('TO',ix_)
        self.xnames=arrset(self.xnames,iv,'TO')
        [iv,ix_,_] = s2mpj_ii('TOB',ix_)
        self.xnames=arrset(self.xnames,iv,'TOB')
        [iv,ix_,_] = s2mpj_ii('NU',ix_)
        self.xnames=arrset(self.xnames,iv,'NU')
        [iv,ix_,_] = s2mpj_ii('D',ix_)
        self.xnames=arrset(self.xnames,iv,'D')
        [iv,ix_,_] = s2mpj_ii('P',ix_)
        self.xnames=arrset(self.xnames,iv,'P')
        [iv,ix_,_] = s2mpj_ii('E',ix_)
        self.xnames=arrset(self.xnames,iv,'E')
        [iv,ix_,_] = s2mpj_ii('P0',ix_)
        self.xnames=arrset(self.xnames,iv,'P0')
        [iv,ix_,_] = s2mpj_ii('G',ix_)
        self.xnames=arrset(self.xnames,iv,'G')
        [iv,ix_,_] = s2mpj_ii('DM',ix_)
        self.xnames=arrset(self.xnames,iv,'DM')
        [iv,ix_,_] = s2mpj_ii('FR',ix_)
        self.xnames=arrset(self.xnames,iv,'FR')
        [iv,ix_,_] = s2mpj_ii('TOLIM',ix_)
        self.xnames=arrset(self.xnames,iv,'TOLIM')
        [iv,ix_,_] = s2mpj_ii('TOBLIM',ix_)
        self.xnames=arrset(self.xnames,iv,'TOBLIM')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('E1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E1')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['F']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('E2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E2')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['K']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('E3',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E3')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['DE']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['D']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['DM']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('E4',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E4')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['DI']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['D']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['DM']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('E5',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E5')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['D']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['P']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['E']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('E6',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E6')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['NU']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['N']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('E7',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E7')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['D']])
        valA = np.append(valA,float(1.5))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['L0']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('E8',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E8')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['L']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['LB']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['FR']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('E9',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E9')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['LB']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('E10',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E10')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['L']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['L0']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['F']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('E11',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E11')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['TO']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('E12',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E12')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['TOB']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('E13',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'E13')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['TO']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['TOLIM']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('E14',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'E14')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['TOB']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['TOBLIM']])
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
        self.gconst = arrset(self.gconst,ig_['E6'],float(-2.0))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        self.xupper[ix_['L0']] = 100.0
        self.xupper[ix_['N']] = 100.0
        self.xupper[ix_['F']] = 30.0
        self.xupper[ix_['K']] = 100.0
        self.xupper[ix_['LB']] = 50.0
        self.xupper[ix_['L']] = 50.0
        self.xupper[ix_['DE']] = 30.0
        self.xupper[ix_['DI']] = 30.0
        self.xupper[ix_['TO']] = 800.0
        self.xupper[ix_['TOB']] = 800.0
        self.xupper[ix_['NU']] = 50.0
        self.xlower[ix_['NU']] = 0.5
        self.xupper[ix_['D']] = 10.0
        self.xlower[ix_['D']] = 0.1
        self.xupper[ix_['P']] = 20.0
        self.xupper[ix_['E']] = 10.0
        self.xupper[ix_['P0']] = 1000.0
        self.xlower[ix_['P0']] = 1.0
        self.xupper[ix_['G']] = 80000.0
        self.xlower[ix_['G']] = 40000.0
        self.xupper[ix_['DM']] = 30.0
        self.xlower[ix_['DM']] = 0.1
        self.xupper[ix_['FR']] = 50.0
        self.xupper[ix_['TOLIM']] = 1000.0
        self.xlower[ix_['TOLIM']] = 100.0
        self.xupper[ix_['TOBLIM']] = 1000.0
        self.xlower[ix_['TOBLIM']] = 100.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        self.x0[ix_['L0']] = float(1.5000e-01)
        self.x0[ix_['N']] = float(2.4079e+01)
        self.x0[ix_['F']] = float(9.2459e-15)
        self.x0[ix_['K']] = float(0.0000)
        self.x0[ix_['LB']] = float(0.0000)
        self.x0[ix_['L']] = float(1.5000e-01)
        self.x0[ix_['DE']] = float(6.8120)
        self.x0[ix_['DI']] = float(6.6120)
        self.x0[ix_['TO']] = float(0.0000)
        self.x0[ix_['TOB']] = float(0.0000)
        self.x0[ix_['NU']] = float(2.2079e+01)
        self.x0[ix_['D']] = float(1.0000e-01)
        self.x0[ix_['P']] = float(6.5268e-01)
        self.x0[ix_['E']] = float(5.5268e-01)
        self.x0[ix_['P0']] = float(6.5887e+02)
        self.x0[ix_['G']] = float(6.5887e+04)
        self.x0[ix_['DM']] = float(6.7120)
        self.x0[ix_['FR']] = float(1.5000e-01)
        self.x0[ix_['TOLIM']] = float(1.0000e+02)
        self.x0[ix_['TOBLIM']] = float(1.0000e+02)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        [it,iet_,_] = s2mpj_ii( 'en311d14', iet_)
        elftv = loaset(elftv,it,0,'V')
        elftv = loaset(elftv,it,1,'W')
        elftv = loaset(elftv,it,2,'X')
        elftv = loaset(elftv,it,3,'Y')
        elftv = loaset(elftv,it,4,'Z')
        [it,iet_,_] = s2mpj_ii( 'en14d31', iet_)
        elftv = loaset(elftv,it,0,'W')
        elftv = loaset(elftv,it,1,'X')
        elftv = loaset(elftv,it,2,'Y')
        elftv = loaset(elftv,it,3,'Z')
        [it,iet_,_] = s2mpj_ii( 'en11d3', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        elftv = loaset(elftv,it,2,'Z')
        [it,iet_,_] = s2mpj_ii( 'en111d2', iet_)
        elftv = loaset(elftv,it,0,'W')
        elftv = loaset(elftv,it,1,'X')
        elftv = loaset(elftv,it,2,'Y')
        elftv = loaset(elftv,it,3,'Z')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        ename = 'EL1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en311d14')
        ielftype = arrset(ielftype,ie,iet_["en311d14"])
        vname = 'DM'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='V')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'NU'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='W')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'P0'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'G'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'D'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='Z')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'EL2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en14d31')
        ielftype = arrset(ielftype,ie,iet_["en14d31"])
        vname = 'G'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='W')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'D'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'DM'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'NU'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='Z')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'EL3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PR')
        ielftype = arrset(ielftype,ie,iet_["en2PR"])
        vname = 'NU'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'P'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'EL4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en11d3')
        ielftype = arrset(ielftype,ie,iet_["en11d3"])
        vname = 'P0'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'DM'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'D'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='Z')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'EL5'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en111d2')
        ielftype = arrset(ielftype,ie,iet_["en111d2"])
        vname = 'G'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='W')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'D'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'E'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'DM'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='Z')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.cupper[np.arange(self.nle)] = np.zeros((self.nle,1))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons   = np.arange(len(self.congrps))
        self.pbclass   = "C-CNLR2-MN-20-14"
        self.objderlvl = 2
        self.conderlvl = [2]


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def e_globs(self):

        import numpy as np
        self.efpar = np.array([]);
        self.efpar = arrset( self.efpar,0,3.1415926535)
        return pbm

    @staticmethod
    def en2PR(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0,0]*EV_[1,0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1,0]
            g_[1] = EV_[0,0]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = 1.0
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def en311d14(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        V3WX = EV_[0,0]**3*EV_[1,0]*EV_[2,0]
        YZ4 = EV_[3,0]*EV_[4,0]**4
        f_   = V3WX/YZ4
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = (3.0*EV_[0,0]**2*EV_[1,0]*EV_[2,0])/YZ4
            g_[1] = (EV_[0,0]**3*EV_[2,0])/YZ4
            g_[2] = (EV_[0,0]**3*EV_[1,0])/YZ4
            g_[3] = -V3WX/(EV_[3,0]*YZ4)
            g_[4] = -(4.0*V3WX)/(EV_[4,0]*YZ4)
            if nargout>2:
                H_ = np.zeros((5,5))
                H_[0,0] = (6.0*EV_[0,0]*EV_[1,0]*EV_[2,0])/YZ4
                H_[0,1] = (3.0*EV_[0,0]**2*EV_[2,0])/YZ4
                H_[1,0] = H_[0,1]
                H_[0,2] = (3.0*EV_[0,0]**2*EV_[1,0])/YZ4
                H_[2,0] = H_[0,2]
                H_[0,3] = -(3.0*EV_[0,0]**2*EV_[1,0])/(YZ4*EV_[3,0])
                H_[3,0] = H_[0,3]
                H_[0,4] = -(12.0*EV_[0,0]**2*EV_[1,0]*EV_[2,0])/(YZ4*EV_[4,0])
                H_[4,0] = H_[0,4]
                H_[1,2] = EV_[0,0]**3/YZ4
                H_[2,1] = H_[1,2]
                H_[1,3] = -(EV_[0,0]**3*EV_[2,0])/(YZ4*EV_[3,0])
                H_[3,1] = H_[1,3]
                H_[1,4] = -(4.0*EV_[0,0]**3*EV_[2,0])/(YZ4*EV_[4,0])
                H_[4,1] = H_[1,4]
                H_[2,3] = -(EV_[0,0]**3*EV_[1,0])/(YZ4*EV_[3,0])
                H_[3,2] = H_[2,3]
                H_[2,4] = -(4.0*EV_[0,0]**3*EV_[1,0])/(YZ4*EV_[4,0])
                H_[4,2] = H_[2,4]
                H_[3,3] = -(2.0*V3WX)/(EV_[3,0]**2*YZ4)
                H_[3,4] = (4.0*V3WX)/(EV_[3,0]*EV_[4,0]*YZ4)
                H_[4,3] = H_[3,4]
                H_[4,4] = (20.0*V3WX)/(EV_[4,0]**2*YZ4)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def en14d31(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        WX4 = EV_[0,0]*EV_[1,0]**4
        Y3Z = EV_[2,0]**3*EV_[3,0]
        f_   = WX4/Y3Z
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1,0]**4/Y3Z
            g_[1] = (4.0*EV_[0,0]*EV_[1,0]**3)/Y3Z
            g_[2] = -(3.0*WX4)/(EV_[2,0]*Y3Z)
            g_[3] = -WX4/(EV_[3,0]*Y3Z)
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,1] = (4.0*EV_[1,0]**3)/Y3Z
                H_[1,0] = H_[0,1]
                H_[0,2] = -(3.0*EV_[1,0]**4)/(EV_[2,0]*Y3Z)
                H_[2,0] = H_[0,2]
                H_[0,3] = -EV_[1,0]**4/(EV_[3,0]*Y3Z)
                H_[3,0] = H_[0,3]
                H_[1,1] = (12.0*EV_[0,0]*EV_[1,0]**2)/Y3Z
                H_[1,2] = -(12.0*EV_[0,0]*EV_[1,0]**3)/(EV_[2,0]*Y3Z)
                H_[2,1] = H_[1,2]
                H_[1,3] = -(4.0*EV_[0,0]*EV_[1,0]**3)/(EV_[3,0]*Y3Z)
                H_[3,1] = H_[1,3]
                H_[2,2] = (12.0*WX4)/(EV_[2,0]**2*Y3Z)
                H_[2,3] = (3.0*WX4)/(EV_[3,0]*EV_[2,0]*Y3Z)
                H_[3,2] = H_[2,3]
                H_[3,3] = (2.0*WX4)/(EV_[3,0]**2*Y3Z)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def en11d3(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = (EV_[0,0]*EV_[1,0])/(self.efpar[0]*EV_[2,0]**3)
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1,0]/(self.efpar[0]*EV_[2,0]**3)
            g_[1] = EV_[0,0]/(self.efpar[0]*EV_[2,0]**3)
            g_[2] = -(3.0*EV_[0,0]*EV_[1,0])/(self.efpar[0]*EV_[2,0]**4)
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = 1.0/(self.efpar[0]*EV_[2,0]**3)
                H_[1,0] = H_[0,1]
                H_[0,2] = -(3.0*EV_[1,0])/(self.efpar[0]*EV_[2,0]**4)
                H_[2,0] = H_[0,2]
                H_[1,2] = -(3.0*EV_[0,0])/(self.efpar[0]*EV_[2,0]**4)
                H_[2,1] = H_[1,2]
                H_[2,2] = (12.0*EV_[0,0]*EV_[1,0])/(self.efpar[0]*EV_[2,0]**5)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def en111d2(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = (EV_[0,0]*EV_[1,0]*EV_[2,0])/(self.efpar[0]*EV_[3,0]**2)
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = (EV_[1,0]*EV_[2,0])/(self.efpar[0]*EV_[3,0]**2)
            g_[1] = (EV_[0,0]*EV_[2,0])/(self.efpar[0]*EV_[3,0]**2)
            g_[2] = (EV_[0,0]*EV_[1,0])/(self.efpar[0]*EV_[3,0]**2)
            g_[3] = -(2.0*EV_[0,0]*EV_[1,0]*EV_[2,0])/(self.efpar[0]*EV_[3,0]**3)
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,1] = EV_[2,0]/(self.efpar[0]*EV_[3,0]**2)
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1,0]/(self.efpar[0]*EV_[3,0]**2)
                H_[2,0] = H_[0,2]
                H_[0,3] = -(2.0*EV_[1,0]*EV_[2,0])/(self.efpar[0]*EV_[3,0]**3)
                H_[3,0] = H_[0,3]
                H_[1,2] = EV_[0,0]/(self.efpar[0]*EV_[3,0]**2)
                H_[2,1] = H_[1,2]
                H_[1,3] = -(2.0*EV_[0,0]*EV_[2,0])/(self.efpar[0]*EV_[3,0]**3)
                H_[3,1] = H_[1,3]
                H_[2,3] = -(2.0*EV_[0,0]*EV_[1,0])/(self.efpar[0]*EV_[3,0]**3)
                H_[3,2] = H_[2,3]
                H_[3,3] = (6.0*EV_[0,0]*EV_[1,0]*EV_[2,0])/(self.efpar[0]*EV_[3,0]**4)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

