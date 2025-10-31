from s2mpjlib import *
class  LUKVLI6(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : LUKVLI6
#    *********
# 
#    Source: Problem 5.6, Generalized Broyden banded function with 
#    exponential constraints, due to L. Luksan and J. Vlcek,
#    "Sparse and partially separable test problems for 
#    unconstrained and equality constrained optimization",
#    Technical Report 767, Inst. Computer Science, Academy of Sciences
#    of the Czech Republic, 182 07 Prague, Czech Republic, 1999
# 
#    Equality constraints changed to inequalities
# 
#    SIF input: Nick Gould, April 2001
# 
#    classification = "C-COOR2-AY-V-V"
# 
#    some useful parameters, including N, the number of variables.
# 
#           Alternative values for the SIF file parameters:
# IE N                   99             $-PARAMETER
# IE N                   999            $-PARAMETER
# IE N                   9999           $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'LUKVLI6'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(9);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
# IE N                   99999          $-PARAMETER
        v_['0'] = 0
        v_['1'] = 1
        v_['2'] = 2
        v_['3'] = 3
        v_['4'] = 4
        v_['5'] = 5
        v_['6'] = 6
        v_['N/2'] = int(np.fix(v_['N']/v_['2']))
        v_['N+1'] = 1+v_['N']
        v_['N-4'] = -4+v_['N']
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
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['I+1'] = 1+I
            v_['I-5'] = -5+I
            [ig,ig_,_] = s2mpj_ii('OBJ'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(I)]])
            valA = np.append(valA,float(2.0))
            v_['A'] = v_['I-5']
            v_['B'] = v_['1']
            v_['A'] = float(v_['A'])
            v_['ABSA'] = np.absolute(v_['A'])
            v_['ABSA'] = int(np.fix(v_['ABSA']))
            v_['B'] = float(v_['B'])
            v_['ABSB'] = np.absolute(v_['B'])
            v_['ABSB'] = int(np.fix(v_['ABSB']))
            v_['ABSA+B'] = v_['ABSA']+v_['ABSB']
            v_['A'] = v_['A']+v_['ABSA+B']
            v_['B'] = v_['B']+v_['ABSA+B']
            v_['A/B'] = int(np.fix(v_['A']/v_['B']))
            v_['B/A'] = int(np.fix(v_['B']/v_['A']))
            v_['SUM'] = v_['A/B']+v_['B/A']
            v_['A'] = v_['A']*v_['A/B']
            v_['B'] = v_['B']*v_['B/A']
            v_['MAXA,B'] = v_['A']+v_['B']
            v_['MAXA,B'] = int(np.fix(v_['MAXA,B']/v_['SUM']))
            v_['MAXA,B'] = v_['MAXA,B']-v_['ABSA+B']
            v_['MAXI-5,1'] = v_['MAXA,B']
            v_['A'] = v_['I+1']
            v_['B'] = v_['N']
            v_['A'] = -1*v_['A']
            v_['B'] = -1*v_['B']
            v_['A'] = float(v_['A'])
            v_['ABSA'] = np.absolute(v_['A'])
            v_['ABSA'] = int(np.fix(v_['ABSA']))
            v_['B'] = float(v_['B'])
            v_['ABSB'] = np.absolute(v_['B'])
            v_['ABSB'] = int(np.fix(v_['ABSB']))
            v_['ABSA+B'] = v_['ABSA']+v_['ABSB']
            v_['A'] = v_['A']+v_['ABSA+B']
            v_['B'] = v_['B']+v_['ABSA+B']
            v_['A/B'] = int(np.fix(v_['A']/v_['B']))
            v_['B/A'] = int(np.fix(v_['B']/v_['A']))
            v_['SUM'] = v_['A/B']+v_['B/A']
            v_['A'] = v_['A']*v_['A/B']
            v_['B'] = v_['B']*v_['B/A']
            v_['MAXA,B'] = v_['A']+v_['B']
            v_['MAXA,B'] = int(np.fix(v_['MAXA,B']/v_['SUM']))
            v_['MINA,B'] = v_['ABSA+B']-v_['MAXA,B']
            v_['MINI+1,N'] = v_['MINA,B']
            for J in range(int(v_['MAXI-5,1']),int(v_['MINI+1,N'])+1):
                [ig,ig_,_] = s2mpj_ii('OBJ'+str(I),ig_)
                gtype = arrset(gtype,ig,'<>')
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(J)]])
                valA = np.append(valA,float(1.0))
        for K in range(int(v_['1']),int(v_['N/2'])+1):
            v_['2K'] = 2*K
            [ig,ig_,_] = s2mpj_ii('C'+str(K),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'C'+str(K))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['2K']))]])
            valA = np.append(valA,float(4.0))
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
        for I in range(int(v_['1']),int(v_['N'])+1):
            self.gconst = arrset(self.gconst,ig_['OBJ'+str(I)],float(-1.0))
        for K in range(int(v_['1']),int(v_['N/2'])+1):
            self.gconst = arrset(self.gconst,ig_['C'+str(K)],float(3.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        for I in range(int(v_['1']),int(v_['N'])+1):
            self.x0[ix_['X'+str(I)]] = float(3.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQR', iet_)
        elftv = loaset(elftv,it,0,'V')
        [it,iet_,_] = s2mpj_ii( 'eCUBE', iet_)
        elftv = loaset(elftv,it,0,'V')
        [it,iet_,_] = s2mpj_ii( 'eXEXP', iet_)
        elftv = loaset(elftv,it,0,'VM')
        elftv = loaset(elftv,it,1,'VP')
        elftv = loaset(elftv,it,2,'V')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['1']),int(v_['N'])+1):
            ename = 'S'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_["eSQR"])
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'C'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eCUBE')
            ielftype = arrset(ielftype,ie,iet_["eCUBE"])
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        for K in range(int(v_['1']),int(v_['N/2'])+1):
            v_['2K'] = 2*K
            v_['2K-1'] = -1+v_['2K']
            v_['2K+1'] = 1+v_['2K']
            ename = 'P'+str(K)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eXEXP')
            ielftype = arrset(ielftype,ie,iet_["eXEXP"])
            vname = 'X'+str(int(v_['2K-1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='VM')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['2K']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['2K+1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='VP')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gL7d3',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['I+1'] = 1+I
            v_['I-5'] = -5+I
            ig = ig_['OBJ'+str(I)]
            self.grftype = arrset(self.grftype,ig,'gL7d3')
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['C'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(5.0))
            v_['A'] = v_['I-5']
            v_['B'] = v_['1']
            v_['A'] = float(v_['A'])
            v_['ABSA'] = np.absolute(v_['A'])
            v_['ABSA'] = int(np.fix(v_['ABSA']))
            v_['B'] = float(v_['B'])
            v_['ABSB'] = np.absolute(v_['B'])
            v_['ABSB'] = int(np.fix(v_['ABSB']))
            v_['ABSA+B'] = v_['ABSA']+v_['ABSB']
            v_['A'] = v_['A']+v_['ABSA+B']
            v_['B'] = v_['B']+v_['ABSA+B']
            v_['A/B'] = int(np.fix(v_['A']/v_['B']))
            v_['B/A'] = int(np.fix(v_['B']/v_['A']))
            v_['SUM'] = v_['A/B']+v_['B/A']
            v_['A'] = v_['A']*v_['A/B']
            v_['B'] = v_['B']*v_['B/A']
            v_['MAXA,B'] = v_['A']+v_['B']
            v_['MAXA,B'] = int(np.fix(v_['MAXA,B']/v_['SUM']))
            v_['MAXA,B'] = v_['MAXA,B']-v_['ABSA+B']
            v_['MAXI-5,1'] = v_['MAXA,B']
            v_['A'] = v_['I+1']
            v_['B'] = v_['N']
            v_['A'] = -1*v_['A']
            v_['B'] = -1*v_['B']
            v_['A'] = float(v_['A'])
            v_['ABSA'] = np.absolute(v_['A'])
            v_['ABSA'] = int(np.fix(v_['ABSA']))
            v_['B'] = float(v_['B'])
            v_['ABSB'] = np.absolute(v_['B'])
            v_['ABSB'] = int(np.fix(v_['ABSB']))
            v_['ABSA+B'] = v_['ABSA']+v_['ABSB']
            v_['A'] = v_['A']+v_['ABSA+B']
            v_['B'] = v_['B']+v_['ABSA+B']
            v_['A/B'] = int(np.fix(v_['A']/v_['B']))
            v_['B/A'] = int(np.fix(v_['B']/v_['A']))
            v_['SUM'] = v_['A/B']+v_['B/A']
            v_['A'] = v_['A']*v_['A/B']
            v_['B'] = v_['B']*v_['B/A']
            v_['MAXA,B'] = v_['A']+v_['B']
            v_['MAXA,B'] = int(np.fix(v_['MAXA,B']/v_['SUM']))
            v_['MINA,B'] = v_['ABSA+B']-v_['MAXA,B']
            v_['MINI+1,N'] = v_['MINA,B']
            for J in range(int(v_['MAXI-5,1']),int(v_['MINI+1,N'])+1):
                ig = ig_['OBJ'+str(I)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['S'+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        for K in range(int(v_['1']),int(v_['N/2'])+1):
            ig = ig_['C'+str(K)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['P'+str(K)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 0.0
#    Solution
# LO SOLTN               6.26382E+04
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle+self.neq,self.m)] = np.zeros((self.nge,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-COOR2-AY-V-V"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSQR(self, nargout,*args):

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
            g_[0] = 2.0*EV_[0]
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
    def eCUBE(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]**3
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 3.0*EV_[0]**2
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 6.0*EV_[0]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eXEXP(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((2,3))
        IV_ = np.zeros(2)
        U_[0,0] = U_[0,0]+1
        U_[0,1] = U_[0,1]-1
        U_[1,0] = U_[1,0]+1
        U_[1,1] = U_[1,1]-1
        U_[1,2] = U_[1,2]-1
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        EXPW = np.exp(IV_[1])
        UEXPW = IV_[0]*EXPW
        f_   = UEXPW
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EXPW
            g_[1] = UEXPW
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = EXPW
                H_[1,0] = H_[0,1]
                H_[1,1] = UEXPW
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gL7d3(self,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        Z = np.absolute(GVAR_)
        f_= Z**(7.0/3.0)
        if nargout>1:
            g_ = 7.0*np.sign(GVAR_)*(Z**(4.0/3.0))/3.0
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 28.0*Z**(1.0/3.0)/9.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

