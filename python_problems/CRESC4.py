from s2mpjlib import *
class  CRESC4(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CRESC4
#    *********
# 
#    This problem consists in finding the crescent of smallest area containing
#    a set of points given in the plane.   This problem arises as a subproblem
#    in pattern recognition and has been suggested by J.P. Rasson.  It
#    originates in the detection of "salt domes" (with the potential presence of
#    oil pockets!) from geological data.
# 
#    The present problem is a simplified version where the crescent is entirely
#    determined by the only four data points.
# 
#    The problem is not convex.
# 
#    A crescent is defined as follows.  Assume one has two circles of center
#    C1 and C2 and of radii r1 and r2 respectively. Assume furthermore that
#    r1 >= r2 and that C2 is within C1.  Assume finally that the distance from
#    C1 to C2 is >= r1 - r2.  Then the crescent is the part of the plane
#    contained in circle 2 but not in circle 1.
# 
#    In order to preserve feasibility at all stages (ensuring that the
#    crescent exists and that its area can be computed), the following
#    parametrization is used:
# 
#    ( C2x, C2y ) = ( C1x, C1y ) + a * d * ( cos(t), sin(t) )
# 
#    r1 = a * d + r 
# 
#    r2 = ( a + 1 ) * d + r
# 
#    with the bounds
# 
#    a >= 1, 0 <= t <= 2 * pi, r2 >= 0 , 0 <= d <= 1.
# 
#    SIF input: Ph. Toint, June 1993.
# 
#    classification = "C-COOR2-MY-6-8"
# 
#    number of points to be included in the crescent.
#    the number of constraints is 2*NP
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'CRESC4'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['NP'] = 4
        v_['X1'] = 1.0
        v_['Y1'] = 0.0
        v_['X2'] = 0.0
        v_['Y2'] = 1.0
        v_['X3'] = 0.0
        v_['Y3'] = -1.0
        v_['X4'] = 0.5
        v_['Y4'] = 0.0
        v_['1'] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [iv,ix_,_] = s2mpj_ii('V1',ix_)
        self.xnames=arrset(self.xnames,iv,'V1')
        [iv,ix_,_] = s2mpj_ii('W1',ix_)
        self.xnames=arrset(self.xnames,iv,'W1')
        [iv,ix_,_] = s2mpj_ii('D',ix_)
        self.xnames=arrset(self.xnames,iv,'D')
        [iv,ix_,_] = s2mpj_ii('A',ix_)
        self.xnames=arrset(self.xnames,iv,'A')
        [iv,ix_,_] = s2mpj_ii('T',ix_)
        self.xnames=arrset(self.xnames,iv,'T')
        [iv,ix_,_] = s2mpj_ii('R',ix_)
        self.xnames=arrset(self.xnames,iv,'R')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['1']),int(v_['NP'])+1):
            [ig,ig_,_] = s2mpj_ii('IS2'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'IS2'+str(I))
            [ig,ig_,_] = s2mpj_ii('OS1'+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'OS1'+str(I))
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
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        self.xlower[ix_['V1']] = -float('Inf')
        self.xupper[ix_['V1']] = +float('Inf')
        self.xlower[ix_['W1']] = -float('Inf')
        self.xupper[ix_['W1']] = +float('Inf')
        self.xlower[ix_['R']] = 0.39
        self.xlower[ix_['A']] = 1.0
        self.xlower[ix_['T']] = 0.0
        self.xupper[ix_['T']] = 6.2831852
        self.xlower[ix_['D']] = 1.0e-8
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        self.x0[ix_['V1']] = float(-40.0)
        self.x0[ix_['W1']] = float(5.0)
        self.x0[ix_['R']] = float(0.75)
        self.x0[ix_['A']] = float(2.0)
        self.x0[ix_['T']] = float(1.5)
        self.x0[ix_['D']] = float(1.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQR1', iet_)
        elftv = loaset(elftv,it,0,'A')
        elftv = loaset(elftv,it,1,'R')
        elftv = loaset(elftv,it,2,'D')
        [it,iet_,_] = s2mpj_ii( 'eSQR2', iet_)
        elftv = loaset(elftv,it,0,'D')
        elftv = loaset(elftv,it,1,'R')
        [it,iet_,_] = s2mpj_ii( 'eSC', iet_)
        elftv = loaset(elftv,it,0,'AZ')
        elftv = loaset(elftv,it,1,'BZ')
        elftv = loaset(elftv,it,2,'DZ')
        [it,iet_,_] = s2mpj_ii( 'eDIST', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftp = []
        elftp = loaset(elftp,it,0,'P')
        [it,iet_,_] = s2mpj_ii( 'eDISTX', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'T')
        elftv = loaset(elftv,it,2,'A')
        elftv = loaset(elftv,it,3,'D')
        elftp = loaset(elftp,it,0,'P')
        [it,iet_,_] = s2mpj_ii( 'eDISTY', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'T')
        elftv = loaset(elftv,it,2,'A')
        elftv = loaset(elftv,it,3,'D')
        elftp = loaset(elftp,it,0,'P')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        ename = 'OB'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSC')
        ielftype = arrset(ielftype,ie,iet_["eSC"])
        vname = 'A'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='AZ')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'R'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='BZ')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'D'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='DZ')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'R2SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQR2')
        ielftype = arrset(ielftype,ie,iet_["eSQR2"])
        vname = 'D'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='D')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'R'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='R')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'R1SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQR1')
        ielftype = arrset(ielftype,ie,iet_["eSQR1"])
        vname = 'D'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='D')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='A')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'R'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='R')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        for I in range(int(v_['1']),int(v_['NP'])+1):
            ename = 'XV1'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eDIST')
            ielftype = arrset(ielftype,ie,iet_["eDIST"])
            vname = 'V1'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='P')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['X'+str(I)]))
            ename = 'XV2'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eDISTX')
            ielftype = arrset(ielftype,ie,iet_["eDISTX"])
            vname = 'V1'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'A'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='A')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'D'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='D')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'T'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='T')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='P')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['X'+str(I)]))
            ename = 'YW1'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eDIST')
            ielftype = arrset(ielftype,ie,iet_["eDIST"])
            vname = 'W1'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='P')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['Y'+str(I)]))
            ename = 'YW2'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eDISTY')
            ielftype = arrset(ielftype,ie,iet_["eDISTY"])
            vname = 'W1'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'A'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='A')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'D'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='D')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'T'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='T')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='P')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['Y'+str(I)]))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        ig = ig_['OBJ']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['OB'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        for I in range(int(v_['1']),int(v_['NP'])+1):
            ig = ig_['IS2'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['XV2'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = posel+1
            self.grelt = loaset(self.grelt,ig,posel,ie_['YW2'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel, 1.)
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['R2SQ'])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
            ig = ig_['OS1'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['XV1'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = posel+1
            self.grelt = loaset(self.grelt,ig,posel,ie_['YW1'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel, 1.)
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['R1SQ'])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution             0.87189692
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.cupper[np.arange(self.nle)] = np.zeros((self.nle,1))
        self.clower[np.arange(self.nle+self.neq,self.m)] = np.zeros((self.nge,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-COOR2-MY-6-8"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSQR1(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        Q = EV_[0]*EV_[2]+EV_[1]
        f_   = Q*Q
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*Q*EV_[2]
            g_[2] = 2.0*Q*EV_[0]
            g_[1] = 2.0*Q
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = 2.0*EV_[2]*EV_[2]
                H_[0,2] = 2.0*(EV_[0]*EV_[2]+Q)
                H_[2,0] = H_[0,2]
                H_[0,1] = 2.0*EV_[2]
                H_[1,0] = H_[0,1]
                H_[2,2] = 2.0*EV_[0]*EV_[0]
                H_[2,1] = 2.0*EV_[0]
                H_[1,2] = H_[2,1]
                H_[1,1] = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eSQR2(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((1,2))
        IV_ = np.zeros(1)
        U_[0,0] = U_[0,0]+1
        U_[0,1] = U_[0,1]+1
        IV_[0] = U_[0:1,:].dot(EV_)
        f_   = IV_[0]*IV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
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

    @staticmethod
    def eSC(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        A = EV_[2]+EV_[1]
        AB = 1.0
        AD = 1.0
        B = EV_[0]*EV_[2]+EV_[1]
        BA = EV_[2]
        BB = 1.0
        BD = EV_[0]
        BAD = 1.0
        D = EV_[0]*EV_[2]
        DA = EV_[2]
        DD = EV_[0]
        DAD = 1.0
        E = 2.0*A*D
        EA = 2.0*A*DA
        EB = 2.0*AB*D
        ED = 2.0*(A*DD+AD*D)
        EAB = 2.0*AB*DA
        EAD = 2.0*(AD*DA+A*DAD)
        EBD = 2.0*AB*DD
        EDD = 2.0*(AD*DD+AD*DD)
        P = 2.0*B*D
        PA = 2.0*(B*DA+BA*D)
        loc_PB = 2.0*BB*D
        PD = 2.0*(B*DD+BD*D)
        PAA = 2.0*(BA*DA+BA*DA)
        PAB = 2.0*BB*DA
        PAD = 2.0*(BD*DA+B*DAD+BAD*D+BA*DD)
        PBD = 2.0*BB*DD
        PDD = 2.0*(BD*DD+BD*DD)
        F = D*D+B*B-A*A
        FA = 2.0*(D*DA+B*BA)
        FB = 2.0*(B*BB-A*AB)
        FD = 2.0*(D*DD+B*BD-A*AD)
        FAA = 2.0*(DA*DA+BA*BA)
        FAB = 2.0*BB*BA
        FAD = 2.0*(DD*DA+D*DAD+BD*BA+B*BAD)
        FBB = 2.0*(BB*BB-AB*AB)
        FBD = 2.0*(BD*BB-AD*AB)
        FDD = 2.0*(DD*DD+BD*BD-AD*AD)
        G = D*D-B*B+A*A
        GA = 2.0*(D*DA-B*BA)
        GB = 2.0*(-B*BB+A*AB)
        GD = 2.0*(D*DD-B*BD+A*AD)
        GAA = 2.0*(DA*DA-BA*BA)
        GAB = -2.0*BB*BA
        GAD = 2.0*(DD*DA+D*DAD-BD*BA-B*BAD)
        GBB = 2.0*(-BB*BB+AB*AB)
        GBD = 2.0*(-BD*BB+AD*AB)
        GDD = 2.0*(DD*DD-BD*BD+AD*AD)
        H = F/P
        I = FA-H*PA
        J = FB-H*loc_PB
        K = FD-H*PD
        HA = I/P
        HB = J/P
        HD = K/P
        IA = FAA-HA*PA-H*PAA
        IB = FAB-HB*PA-H*PAB
        ID = FAD-HD*PA-H*PAD
        JB = FBB-HB*loc_PB
        JD = FBD-HD*loc_PB-H*PBD
        KD = FDD-HD*PD-H*PDD
        HAA = (IA-HA*PA)/P
        HAB = (IB-HA*loc_PB)/P
        HAD = (ID-HA*PD)/P
        HBB = (JB-HB*loc_PB)/P
        HBD = (JD-HB*PD)/P
        HDD = (KD-HD*PD)/P
        L = -G/E
        M = -GA-L*EA
        N = -GB-L*EB
        O = -GD-L*ED
        LA = M/E
        LB = N/E
        LD = O/E
        MA = -GAA-LA*EA
        MB = -GAB-LB*EA-L*EAB
        MD = -GAD-LD*EA-L*EAD
        NB = -GBB-LB*EB
        ND = -GBD-LD*EB-L*EBD
        OD = -GDD-LD*ED-L*EDD
        LAA = (MA-LA*EA)/E
        LAB = (MB-LA*EB)/E
        LAD = (MD-LA*ED)/E
        LBB = (NB-LB*EB)/E
        LBD = (ND-LB*ED)/E
        LDD = (OD-LD*ED)/E
        C = np.arccos(H)
        CH = -1.0/np.sqrt(1.0-H*H)
        CHH = CH*H/(1.0-H*H)
        CA = CH*HA
        CB = CH*HB
        CD = CH*HD
        CAA = CHH*HA*HA+CH*HAA
        CAB = CHH*HA*HB+CH*HAB
        CAD = CHH*HA*HD+CH*HAD
        CBB = CHH*HB*HB+CH*HBB
        CBD = CHH*HB*HD+CH*HBD
        CDD = CHH*HD*HD+CH*HDD
        Q = np.arccos(L)
        QL = -1.0/np.sqrt(1.0-L*L)
        QLL = QL*L/(1.0-L*L)
        QA = QL*LA
        QB = QL*LB
        QD = QL*LD
        QAA = QLL*LA*LA+QL*LAA
        QAB = QLL*LA*LB+QL*LAB
        QAD = QLL*LA*LD+QL*LAD
        QBB = QLL*LB*LB+QL*LBB
        QBD = QLL*LB*LD+QL*LBD
        QDD = QLL*LD*LD+QL*LDD
        R = B*B*C
        RA = 2.0*B*BA*C+B*B*CA
        RB = 2.0*B*BB*C+B*B*CB
        RD = 2.0*B*BD*C+B*B*CD
        RAA = 2.0*(BA*BA*C+B*BA*CA+B*BA*CA)+B*B*CAA
        RAB = 2.0*(BB*BA*C+B*BA*CB+B*BB*CA)+B*B*CAB
        RAD = 2.0*(BD*BA*C+B*BAD*C+B*BA*CD+B*BD*CA)+B*B*CAD
        RBB = 2.0*(BB*BB*C+B*BB*CB+B*BB*CB)+B*B*CBB
        RBD = 2.0*(BD*BB*C+B*BB*CD+B*BD*CB)+B*B*CBD
        RDD = 2.0*(BD*BD*C+B*BD*CD+B*BD*CD)+B*B*CDD
        S = A*A*Q
        SA = A*A*QA
        SB = 2.0*A*AB*Q+A*A*QB
        SD = 2.0*A*AD*Q+A*A*QD
        SAA = A*A*QAA
        SAB = 2.0*A*AB*QA+A*A*QAB
        SAD = 2.0*A*AD*QA+A*A*QAD
        SBB = 2.0*(AB*AB*Q+A*AB*QB+A*AB*QB)+A*A*QBB
        SBD = 2.0*(AD*AB*Q+A*AB*QD+A*AD*QB)+A*A*QBD
        SDD = 2.0*(AD*AD*Q+A*AD*QD+A*AD*QD)+A*A*QDD
        SQ = np.sin(Q)
        CQ = L
        W = 0.5*E*SQ
        WA = 0.5*(EA*SQ+E*CQ*QA)
        WB = 0.5*(EB*SQ+E*CQ*QB)
        WD = 0.5*(ED*SQ+E*CQ*QD)
        WAA = 0.5*(EA*CQ*QA+EA*CQ*QA-E*SQ*QA*QA+E*CQ*QAA)
        WAB = 0.5*(EAB*SQ+EA*CQ*QB+EB*CQ*QA-E*SQ*QB*QA+E*CQ*QAB)
        WAD = 0.5*(EAD*SQ+EA*CQ*QD+ED*CQ*QA-E*SQ*QD*QA+E*CQ*QAD)
        WBB = 0.5*(EB*CQ*QB+EB*CQ*QB-E*SQ*QB*QB+E*CQ*QBB)
        WBD = 0.5*(EBD*SQ+EB*CQ*QD+ED*CQ*QB-E*SQ*QD*QB+E*CQ*QBD)
        WDD = 0.5*(EDD*SQ+ED*CQ*QD+ED*CQ*QD-E*SQ*QD*QD+E*CQ*QDD)
        V = S-R+W
        VA = SA-RA+WA
        VB = SB-RB+WB
        VD = SD-RD+WD
        VAA = SAA-RAA+WAA
        VAB = SAB-RAB+WAB
        VAD = SAD-RAD+WAD
        VBB = SBB-RBB+WBB
        VBD = SBD-RBD+WBD
        VDD = SDD-RDD+WDD
        f_   = V
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = VA
            g_[1] = VB
            g_[2] = VD
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = VAA
                H_[0,1] = VAB
                H_[1,0] = H_[0,1]
                H_[0,2] = VAD
                H_[2,0] = H_[0,2]
                H_[1,1] = VBB
                H_[1,2] = VBD
                H_[2,1] = H_[1,2]
                H_[2,2] = VDD
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eDIST(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = (EV_[0]-self.elpar[iel_][0])**2
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*(EV_[0]-self.elpar[iel_][0])
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
    def eDISTY(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        ST = np.sin(EV_[1])
        CT = np.cos(EV_[1])
        B = EV_[0]+EV_[2]*EV_[3]*ST-self.elpar[iel_][0]
        BA = EV_[3]*ST
        BD = EV_[2]*ST
        BT = EV_[2]*EV_[3]*CT
        f_   = B*B
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*B
            g_[2] = 2.0*B*BA
            g_[3] = 2.0*B*BD
            g_[1] = 2.0*B*BT
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,0] = 2.0
                H_[0,2] = 2.0*BA
                H_[2,0] = H_[0,2]
                H_[0,3] = 2.0*BD
                H_[3,0] = H_[0,3]
                H_[0,1] = 2.0*BT
                H_[1,0] = H_[0,1]
                H_[2,2] = 2.0*BA*BA
                H_[2,3] = 2.0*(BD*BA+B*ST)
                H_[3,2] = H_[2,3]
                H_[2,1] = 2.0*(BT*BA+B*EV_[3]*CT)
                H_[1,2] = H_[2,1]
                H_[3,3] = 2.0*BD*BD
                H_[3,1] = 2.0*(BT*BD+B*EV_[2]*CT)
                H_[1,3] = H_[3,1]
                H_[1,1] = 2.0*(BT*BT-B*EV_[2]*EV_[3]*ST)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eDISTX(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        ST = np.sin(EV_[1])
        CT = np.cos(EV_[1])
        B = EV_[0]+EV_[2]*EV_[3]*CT-self.elpar[iel_][0]
        BA = EV_[3]*CT
        BD = EV_[2]*CT
        BT = -EV_[2]*EV_[3]*ST
        f_   = B*B
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*B
            g_[2] = 2.0*B*BA
            g_[3] = 2.0*B*BD
            g_[1] = 2.0*B*BT
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,0] = 2.0
                H_[0,2] = 2.0*BA
                H_[2,0] = H_[0,2]
                H_[0,3] = 2.0*BD
                H_[3,0] = H_[0,3]
                H_[0,1] = 2.0*BT
                H_[1,0] = H_[0,1]
                H_[2,2] = 2.0*BA*BA
                H_[2,3] = 2.0*(BD*BA+B*CT)
                H_[3,2] = H_[2,3]
                H_[2,1] = 2.0*(BT*BA-B*EV_[3]*ST)
                H_[1,2] = H_[2,1]
                H_[3,3] = 2.0*BD*BD
                H_[3,1] = 2.0*(BT*BD-B*EV_[2]*ST)
                H_[1,3] = H_[3,1]
                H_[1,1] = 2.0*(BT*BT-B*EV_[2]*EV_[3]*CT)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

