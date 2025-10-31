from s2mpjlib import *
class  DEVGLA2NE(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DEVGLA2NE
#    *********
# 
#    SCIPY global optimization benchmark example DeVilliersGlasser02
# 
#    Fit: y  = x_1 x_2^t tanh ( t x_3 + sin( t x_4 ) ) cos( t e^x_5 )  +  e
# 
#    Source:  Problem from the SCIPY benchmark set
#      https://github.com/scipy/scipy/tree/master/benchmarks/ ...
#              benchmarks/go_benchmark_functions
# 
#    Nonlinear-equation formulation of DEVGLA2.SIF
# 
#    SIF input: Nick Gould, Jan 2020
# 
#    classification = "C-CNOR2-MN-5-16"
# 
#    Number of data values
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'DEVGLA2NE'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['M'] = 16
        v_['N'] = 5
        v_['1'] = 1
        v_['A'] = 1.27
        v_['LNA'] = np.log(v_['A'])
        for I in range(int(v_['1']),int(v_['M'])+1):
            v_['RI'] = float(I)
            v_['RIM1'] = -1.0+v_['RI']
            v_['T'] = 0.1*v_['RIM1']
            v_['T'+str(I)] = v_['T']
            v_['TLNA'] = v_['T']*v_['LNA']
            v_['AT'] = np.exp(v_['TLNA'])
            v_['TP'] = 3.012*v_['T']
            v_['TP2'] = 2.13*v_['T']
            v_['STP2'] = np.sin(v_['TP2'])
            v_['TPA'] = v_['TP']+v_['STP2']
            v_['HTPA'] = np.tanh(v_['TPA'])
            v_['EC'] = np.exp(0.507)
            v_['ECT'] = v_['EC']*v_['T']
            v_['CECT'] = np.cos(v_['ECT'])
            v_['P'] = v_['AT']*v_['HTPA']
            v_['PP'] = v_['P']*v_['CECT']
            v_['PPP'] = 53.81*v_['PP']
            v_['Y'+str(I)] = v_['PPP']
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
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['1']),int(v_['M'])+1):
            [ig,ig_,_] = s2mpj_ii('F'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'F'+str(I))
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
        for I in range(int(v_['1']),int(v_['M'])+1):
            self.gconst = arrset(self.gconst,ig_['F'+str(I)],float(v_['Y'+str(I)]))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        if('X1' in ix_):
            self.x0[ix_['X1']] = float(20.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X1']),float(20.0)))
        if('X2' in ix_):
            self.x0[ix_['X2']] = float(2.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X2']),float(2.0)))
        if('X3' in ix_):
            self.x0[ix_['X3']] = float(2.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X3']),float(2.0)))
        if('X4' in ix_):
            self.x0[ix_['X4']] = float(2.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X4']),float(2.0)))
        if('X5' in ix_):
            self.x0[ix_['X5']] = float(0.2)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X5']),float(0.2)))
        pass
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eDG2', iet_)
        elftv = loaset(elftv,it,0,'X1')
        elftv = loaset(elftv,it,1,'X2')
        elftv = loaset(elftv,it,2,'X3')
        elftv = loaset(elftv,it,3,'X4')
        elftv = loaset(elftv,it,4,'X5')
        elftp = []
        elftp = loaset(elftp,it,0,'T')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['1']),int(v_['M'])+1):
            ename = 'E'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eDG2')
            ielftype = arrset(ielftype,ie,iet_["eDG2"])
            vname = 'X1'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='X1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='X2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='X3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X4'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='X4')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X5'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='X5')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='T')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['T'+str(I)]))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['M'])+1):
            ig = ig_['F'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        self.objlower = 0.0
#    Solution
# LO SOLUTION            0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CNOR2-MN-5-16"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eDG2(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        X2T = EV_[1]**self.elpar[iel_][0]
        F2 = X2T
        F2X2 = self.elpar[iel_][0]*EV_[1]**(self.elpar[iel_][0]-1.0e0)
        F2X2X2  = (
              self.elpar[iel_][0]*(self.elpar[iel_][0]-1.0e0)*EV_[1]**(self.elpar[iel_][0]-2.0e0))
        X3T = EV_[2]*self.elpar[iel_][0]
        X4T = EV_[3]*self.elpar[iel_][0]
        SINX4T = np.sin(X4T)
        COSX4T = np.cos(X4T)
        A = X3T+SINX4T
        AX3 = self.elpar[iel_][0]
        AX4 = self.elpar[iel_][0]*COSX4T
        AX4X4 = -self.elpar[iel_][0]*self.elpar[iel_][0]*SINX4T
        F = np.tanh(A)
        FA = 1.0/(np.cosh(A))**2
        FAA = -2.0*FA*F
        F3 = F
        F3X3 = FA*AX3
        F3X4 = FA*AX4
        F3X3X3 = FAA*AX3*AX3
        F3X3X4 = FAA*AX3*AX4
        F3X4X4 = FA*AX4X4+FAA*AX4*AX4
        EX5 = np.exp(EV_[4])
        TEX5 = self.elpar[iel_][0]*EX5
        STEX5 = np.sin(TEX5)
        CTEX5 = np.cos(TEX5)
        F4 = CTEX5
        F4X5 = -STEX5*TEX5
        F4X5X5 = -STEX5*TEX5-CTEX5*TEX5*TEX5
        f_   = EV_[0]*F2*F3*F4
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = F2*F3*F4
            g_[1] = EV_[0]*F2X2*F3*F4
            g_[2] = EV_[0]*F2*F3X3*F4
            g_[3] = EV_[0]*F2*F3X4*F4
            g_[4] = EV_[0]*F2*F3*F4X5
            if nargout>2:
                H_ = np.zeros((5,5))
                H_[0,1] = F2X2*F3*F4
                H_[1,0] = H_[0,1]
                H_[0,2] = F2*F3X3*F4
                H_[2,0] = H_[0,2]
                H_[0,3] = F2*F3X4*F4
                H_[3,0] = H_[0,3]
                H_[0,4] = F2*F3*F4X5
                H_[4,0] = H_[0,4]
                H_[1,1] = EV_[0]*F2X2X2*F3*F4
                H_[1,2] = EV_[0]*F2X2*F3X3*F4
                H_[2,1] = H_[1,2]
                H_[1,3] = EV_[0]*F2X2*F3X4*F4
                H_[3,1] = H_[1,3]
                H_[1,4] = EV_[0]*F2X2*F3*F4X5
                H_[4,1] = H_[1,4]
                H_[2,2] = EV_[0]*F2*F3X3X3*F4
                H_[2,3] = EV_[0]*F2*F3X3X4*F4
                H_[3,2] = H_[2,3]
                H_[2,4] = EV_[0]*F2*F3X3*F4X5
                H_[4,2] = H_[2,4]
                H_[3,3] = EV_[0]*F2*F3X4X4*F4
                H_[3,4] = EV_[0]*F2*F3X4*F4X5
                H_[4,3] = H_[3,4]
                H_[4,4] = EV_[0]*F2*F3*F4X5X5
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

