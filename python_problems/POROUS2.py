from s2mpjlib import *
class  POROUS2(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : POROUS2
#    *********
# 
#    The problem is to solve the porous medium equation on the unit square.
#    The equation is
# 
#        \Delta ( u^2 ) + d \frac{\partial}{\partial x_1}( u^3 ) + f = 0
# 
#    within the domain.  The boundary condition are that u = 1 on the bottom
#    and left sides and u = 0 on the top and right sides.  Discretization is
#    using the usual central differences. The function f is a point source of
#    maginitude 50 at the lower left grid point.  The initial approximation
#    is a discretization of 1 - x_1 x_2.
# 
#    Source: example 3.2.4 in
#    S. Eisenstat and H. Walker,
#    "Choosing the forcing terms in an inexact Newton method"
#    Report 6/94/75, Dept of Maths, Utah State University, 1994.
# 
#    SIF input: Ph. Toint, July 1994.
# 
#    classification = "C-CNOR2-MN-V-V"
# 
#    P is the number of points in one side of the unit square.
#    There are P*P variables.
# 
#           Alternative values for the SIF file parameters:
# IE P                   32             $-PARAMETER      original value
# IE P                   64             $-PARAMETER 
# IE P                   72             $-PARAMETER 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'POROUS2'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['P'] = int(5);  #  SIF file default value
        else:
            v_['P'] = int(args[0])
        if nargin<2:
            v_['D'] = float(-50.0);  #  SIF file default value
        else:
            v_['D'] = float(args[1])
        v_['1'] = 1
        v_['2'] = 2
        v_['P-1'] = -1+v_['P']
        v_['RP-1'] = float(v_['P-1'])
        v_['H'] = 1.0/v_['RP-1']
        v_['H2'] = v_['H']*v_['H']
        v_['1/H2'] = 1.0/v_['H2']
        v_['2H'] = 2.0*v_['H']
        v_['D/2H'] = v_['D']/v_['2H']
        v_['-D/2H'] = -1.0*v_['D/2H']
        v_['-4/H2'] = -4.0*v_['1/H2']
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
                [iv,ix_,_] = s2mpj_ii('U'+str(I)+','+str(J),ix_)
                self.xnames=arrset(self.xnames,iv,'U'+str(I)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['2']),int(v_['P-1'])+1):
            for J in range(int(v_['2']),int(v_['P-1'])+1):
                [ig,ig_,_] = s2mpj_ii('G'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'G'+str(I)+','+str(J))
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
        self.gconst  = (
              arrset(self.gconst,ig_['G'+str(int(v_['P-1']))+','+str(int(v_['P-1']))],float(-50.0)))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        for J in range(int(v_['1']),int(v_['P'])+1):
            self.xlower[ix_['U'+str(int(v_['1']))+','+str(J)]] = 1.0
            self.xupper[ix_['U'+str(int(v_['1']))+','+str(J)]] = 1.0
            self.xlower[ix_['U'+str(int(v_['P']))+','+str(J)]] = 0.0
            self.xupper[ix_['U'+str(int(v_['P']))+','+str(J)]] = 0.0
        for I in range(int(v_['2']),int(v_['P-1'])+1):
            self.xlower[ix_['U'+str(I)+','+str(int(v_['P']))]] = 1.0
            self.xupper[ix_['U'+str(I)+','+str(int(v_['P']))]] = 1.0
            self.xlower[ix_['U'+str(I)+','+str(int(v_['1']))]] = 0.0
            self.xupper[ix_['U'+str(I)+','+str(int(v_['1']))]] = 0.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        for I in range(int(v_['1']),int(v_['P'])+1):
            for J in range(int(v_['1']),int(v_['P'])+1):
                v_['RI'] = float(I)
                v_['RJ'] = float(J)
                v_['I-1'] = -1.0+v_['RI']
                v_['J-1'] = -1.0+v_['RJ']
                v_['X1'] = v_['I-1']*v_['H']
                v_['X2'] = v_['J-1']*v_['H']
                v_['X1X2'] = v_['X1']*v_['X2']
                v_['MX1X2'] = -1.0*v_['X1X2']
                v_['UIJ'] = 1.0+v_['MX1X2']
                self.x0[ix_['U'+str(I)+','+str(J)]] = float(v_['UIJ'])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'U')
        [it,iet_,_] = s2mpj_ii( 'eCB', iet_)
        elftv = loaset(elftv,it,0,'U')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['1']),int(v_['P'])+1):
            for J in range(int(v_['1']),int(v_['P'])+1):
                ename = 'US'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eSQ')
                ielftype = arrset(ielftype,ie,iet_["eSQ"])
                vname = 'U'+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='U')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                ename = 'UC'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eCB')
                ielftype = arrset(ielftype,ie,iet_["eCB"])
                vname = 'U'+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='U')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['2']),int(v_['P-1'])+1):
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            for J in range(int(v_['2']),int(v_['P-1'])+1):
                v_['J-1'] = -1+J
                v_['J+1'] = 1+J
                ig = ig_['G'+str(I)+','+str(J)]
                posel = len(self.grelt[ig])
                self.grelt  = (
                      loaset(self.grelt,ig,posel,ie_['US'+str(int(v_['I+1']))+','+str(J)]))
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(v_['1/H2']))
                posel = len(self.grelt[ig])
                self.grelt  = (
                      loaset(self.grelt,ig,posel,ie_['US'+str(int(v_['I-1']))+','+str(J)]))
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(v_['1/H2']))
                posel = len(self.grelt[ig])
                self.grelt  = (
                      loaset(self.grelt,ig,posel,ie_['US'+str(I)+','+str(int(v_['J-1']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(v_['1/H2']))
                posel = len(self.grelt[ig])
                self.grelt  = (
                      loaset(self.grelt,ig,posel,ie_['US'+str(I)+','+str(int(v_['J+1']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(v_['1/H2']))
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['US'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(v_['-4/H2']))
                posel = len(self.grelt[ig])
                self.grelt  = (
                      loaset(self.grelt,ig,posel,ie_['UC'+str(int(v_['I+1']))+','+str(J)]))
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(v_['D/2H']))
                posel = len(self.grelt[ig])
                self.grelt  = (
                      loaset(self.grelt,ig,posel,ie_['UC'+str(int(v_['I-1']))+','+str(J)]))
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(v_['-D/2H']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 0.0
#    Solution
# LO SOLTN               0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CNOR2-MN-V-V"
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

    @staticmethod
    def eCB(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0,0]*EV_[0,0]*EV_[0,0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 3.0*EV_[0,0]*EV_[0,0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 6.0*EV_[0,0]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

