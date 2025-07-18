from s2mpjlib import *
class  ORTHREGF(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : ORTHREGF
#    *********
# 
#    An orthogonal regression problem
# 
#    The problem is to fit (orthogonally) an torus to a
#    set of points in 3D space. This set of points is generated by
#    perturbing a first set lying exactly on a predefined torus
#    centered at the origin.
# 
#    Source:
#    M. Gulliksson,
#    "Algorithms for nonlinear Least-squares with Applications to
#    Orthogonal Regression",
#    UMINF-178.90, University of Umea, Sweden, 1990.
# 
#    SIF input: Ph. Toint, June 1990.
#               minor correction by Ph. Shott, Jan 1995.
# 
#    classification = "C-CQOR2-AY-V-V"
# 
#    square root of the number of data points
#    (number of variables = 3 * NPTS**2 + 5 )
# 
#           Alternative values for the SIF file parameters:
# IE NPTS                5              $-PARAMETER n = 80    original value
# IE NPTS                7              $-PARAMETER n = 152
# IE NPTS                10             $-PARAMETER n = 305
# IE NPTS                15             $-PARAMETER n = 680
# IE NPTS                20             $-PARAMETER n = 1205
# IE NPTS                40             $-PARAMETER n = 4805
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 21 VI 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'ORTHREGF'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['NPTS'] = int(5);  #  SIF file default value
        else:
            v_['NPTS'] = int(args[0])
        v_['TP4'] = 1.7
        v_['TP5'] = 0.8
        v_['PSEED'] = 237.1531
        v_['PSIZE'] = 0.2
        v_['1'] = 1
        v_['5'] = 5
        v_['PI'] = 3.1415926535
        v_['2PI'] = 2.0*v_['PI']
        v_['RNPTS'] = float(v_['NPTS'])
        v_['ICR0'] = 1.0/v_['RNPTS']
        v_['INCR'] = v_['ICR0']*v_['2PI']
        for I in range(int(v_['1']),int(v_['NPTS'])+1):
            v_['I-1'] = -1+I
            v_['RI-1'] = float(v_['I-1'])
            v_['THETA1'] = v_['RI-1']*v_['INCR']
            v_['ST1'] = np.sin(v_['THETA1'])
            v_['CT1'] = np.cos(v_['THETA1'])
            v_['P5CT1'] = v_['TP5']*v_['CT1']
            v_['P4P5CT1'] = v_['TP4']+v_['P5CT1']
            v_['R3'] = v_['TP5']*v_['ST1']
            for J in range(int(v_['1']),int(v_['NPTS'])+1):
                v_['J-1'] = -1+J
                v_['RJ-1'] = float(v_['J-1'])
                v_['THETA2'] = v_['RJ-1']*v_['INCR']
                v_['ST2'] = np.sin(v_['THETA2'])
                v_['CT2'] = np.cos(v_['THETA2'])
                v_['R1'] = v_['P4P5CT1']*v_['CT2']
                v_['R2'] = v_['P4P5CT1']*v_['ST2']
                v_['XSEED'] = v_['THETA2']*v_['PSEED']
                v_['SSEED'] = np.cos(v_['XSEED'])
                v_['PER-1'] = v_['PSIZE']*v_['SSEED']
                v_['PERT'] = 1.0+v_['PER-1']
                v_['XD'+str(I)+','+str(J)] = v_['R1']*v_['PERT']
                v_['YD'+str(I)+','+str(J)] = v_['R2']*v_['PERT']
                v_['ZD'+str(I)+','+str(J)] = v_['R3']*v_['PERT']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['5'])+1):
            [iv,ix_,_] = s2mpj_ii('P'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'P'+str(I))
        for I in range(int(v_['1']),int(v_['NPTS'])+1):
            for J in range(int(v_['1']),int(v_['NPTS'])+1):
                [iv,ix_,_] = s2mpj_ii('X'+str(I)+','+str(J),ix_)
                self.xnames=arrset(self.xnames,iv,'X'+str(I)+','+str(J))
                [iv,ix_,_] = s2mpj_ii('Y'+str(I)+','+str(J),ix_)
                self.xnames=arrset(self.xnames,iv,'Y'+str(I)+','+str(J))
                [iv,ix_,_] = s2mpj_ii('Z'+str(I)+','+str(J),ix_)
                self.xnames=arrset(self.xnames,iv,'Z'+str(I)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['1']),int(v_['NPTS'])+1):
            for J in range(int(v_['1']),int(v_['NPTS'])+1):
                [ig,ig_,_] = s2mpj_ii('OX'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(I)+','+str(J)]])
                valA = np.append(valA,float(1.0))
                [ig,ig_,_] = s2mpj_ii('OY'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['Y'+str(I)+','+str(J)]])
                valA = np.append(valA,float(1.0))
                [ig,ig_,_] = s2mpj_ii('OZ'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['Z'+str(I)+','+str(J)]])
                valA = np.append(valA,float(1.0))
                [ig,ig_,_] = s2mpj_ii('A'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'A'+str(I)+','+str(J))
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
        for I in range(int(v_['1']),int(v_['NPTS'])+1):
            for J in range(int(v_['1']),int(v_['NPTS'])+1):
                self.gconst  = (
                      arrset(self.gconst,ig_['OX'+str(I)+','+str(J)],float(v_['XD'+str(I)+','+str(J)])))
                self.gconst  = (
                      arrset(self.gconst,ig_['OY'+str(I)+','+str(J)],float(v_['YD'+str(I)+','+str(J)])))
                self.gconst  = (
                      arrset(self.gconst,ig_['OZ'+str(I)+','+str(J)],float(v_['ZD'+str(I)+','+str(J)])))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        self.xlower[ix_['P4']] = 0.001
        self.xlower[ix_['P5']] = 0.001
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        if('P1' in ix_):
            self.x0[ix_['P1']] = float(1.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P1']),float(1.0)))
        if('P2' in ix_):
            self.x0[ix_['P2']] = float(0.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P2']),float(0.0)))
        if('P3' in ix_):
            self.x0[ix_['P3']] = float(1.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P3']),float(1.0)))
        if('P4' in ix_):
            self.x0[ix_['P4']] = float(1.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P4']),float(1.0)))
        if('P5' in ix_):
            self.x0[ix_['P5']] = float(0.5)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['P5']),float(0.5)))
        for I in range(int(v_['1']),int(v_['NPTS'])+1):
            for J in range(int(v_['1']),int(v_['NPTS'])+1):
                if('X'+str(I)+','+str(J) in ix_):
                    self.x0[ix_['X'+str(I)+','+str(J)]] = float(v_['XD'+str(I)+','+str(J)])
                else:
                    self.y0  = (
                          arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X'+str(I)+','+str(J)]),float(v_['XD'+str(I)+','+str(J)])))
                if('Y'+str(I)+','+str(J) in ix_):
                    self.x0[ix_['Y'+str(I)+','+str(J)]] = float(v_['YD'+str(I)+','+str(J)])
                else:
                    self.y0  = (
                          arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Y'+str(I)+','+str(J)]),float(v_['YD'+str(I)+','+str(J)])))
                if('Z'+str(I)+','+str(J) in ix_):
                    self.x0[ix_['Z'+str(I)+','+str(J)]] = float(v_['ZD'+str(I)+','+str(J)])
                else:
                    self.y0  = (
                          arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Z'+str(I)+','+str(J)]),float(v_['ZD'+str(I)+','+str(J)])))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eTA', iet_)
        elftv = loaset(elftv,it,0,'XX')
        elftv = loaset(elftv,it,1,'YY')
        elftv = loaset(elftv,it,2,'A')
        elftv = loaset(elftv,it,3,'B')
        elftv = loaset(elftv,it,4,'C')
        [it,iet_,_] = s2mpj_ii( 'eISQ', iet_)
        elftv = loaset(elftv,it,0,'Z')
        elftv = loaset(elftv,it,1,'P')
        [it,iet_,_] = s2mpj_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'XX')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['1']),int(v_['NPTS'])+1):
            for J in range(int(v_['1']),int(v_['NPTS'])+1):
                ename = 'EA'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eTA')
                ielftype = arrset(ielftype,ie,iet_["eTA"])
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                posev = np.where(elftv[ielftype[ie]]=='XX')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'Y'+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                posev = np.where(elftv[ielftype[ie]]=='YY')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'P1'
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                posev = np.where(elftv[ielftype[ie]]=='A')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'P2'
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                posev = np.where(elftv[ielftype[ie]]=='B')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'P4'
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                posev = np.where(elftv[ielftype[ie]]=='C')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                ename = 'EB'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eISQ')
                ielftype = arrset(ielftype,ie,iet_["eISQ"])
                vname = 'Z'+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                posev = np.where(elftv[ielftype[ie]]=='Z')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'P3'
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                posev = np.where(elftv[ielftype[ie]]=='P')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                ename = 'EC'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eSQ')
                ielftype = arrset(ielftype,ie,iet_["eSQ"])
                vname = 'P5'
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                posev = np.where(elftv[ielftype[ie]]=='XX')[0]
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
        for I in range(int(v_['1']),int(v_['NPTS'])+1):
            for J in range(int(v_['1']),int(v_['NPTS'])+1):
                ig = ig_['OX'+str(I)+','+str(J)]
                self.grftype = arrset(self.grftype,ig,'gL2')
                ig = ig_['OY'+str(I)+','+str(J)]
                self.grftype = arrset(self.grftype,ig,'gL2')
                ig = ig_['OZ'+str(I)+','+str(J)]
                self.grftype = arrset(self.grftype,ig,'gL2')
                ig = ig_['A'+str(I)+','+str(J)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['EA'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,1.)
                posel = posel+1
                self.grelt = loaset(self.grelt,ig,posel,ie_['EB'+str(I)+','+str(J)])
                self.grelw = loaset(self.grelw,ig,posel, 1.)
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['EC'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 0.0
#    Solution
# LO SOLTN(5)            0.990089426
# LO SOLTN(7)            1.315031322
# LO SOLTN(10)           4.515848902
# LO SOLTN(15)           9.185538338
# LO SOLTN(20)           16.20054380
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
        self.pbclass   = "C-CQOR2-AY-V-V"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eTA(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((3,5))
        IV_ = np.zeros(3)
        U_[0,0] = U_[0,0]+1
        U_[0,2] = U_[0,2]-1
        U_[1,1] = U_[1,1]+1
        U_[1,3] = U_[1,3]-1
        U_[2,4] = U_[2,4]+1
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        IV_[2] = U_[2:3,:].dot(EV_)
        CCSQ = IV_[2]*IV_[2]
        CCCB = CCSQ*IV_[2]
        XXYY = IV_[0]*IV_[0]+IV_[1]*IV_[1]
        T = XXYY/CCSQ
        DTDX = 2.0*IV_[0]/CCSQ
        DTDY = 2.0*IV_[1]/CCSQ
        DTDC = -2.0*XXYY/CCCB
        D2TDX2 = 2.0/CCSQ
        D2TDY2 = 2.0/CCSQ
        D2TDC2 = 6.0*XXYY/(CCSQ*CCSQ)
        D2TDXC = -4.0*IV_[0]/CCCB
        D2TDYC = -4.0*IV_[1]/CCCB
        S = np.sqrt(T)
        R = 0.5/S
        DSDX = R*DTDX
        DSDY = R*DTDY
        DSDC = R*DTDC
        D2SDX2 = R*(D2TDX2-0.5*DTDX*DTDX/T)
        D2SDY2 = R*(D2TDY2-0.5*DTDY*DTDY/T)
        D2SDC2 = R*(D2TDC2-0.5*DTDC*DTDC/T)
        D2SDXY = -0.5*DTDX*DSDY/T
        D2SDXC = R*(D2TDXC-0.5*DTDX*DTDC/T)
        D2SDYC = R*(D2TDYC-0.5*DTDY*DTDC/T)
        SS = S-1.0
        SPS = SS+SS
        f_   = SS*SS
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = SPS*DSDX
            g_[1] = SPS*DSDY
            g_[2] = SPS*DSDC
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = SPS*D2SDX2+2.0*DSDX*DSDX
                H_[0,1] = SPS*D2SDXY+2.0*DSDX*DSDY
                H_[1,0] = H_[0,1]
                H_[0,2] = SPS*D2SDXC+2.0*DSDX*DSDC
                H_[2,0] = H_[0,2]
                H_[1,1] = SPS*D2SDY2+2.0*DSDY*DSDY
                H_[1,2] = SPS*D2SDYC+2.0*DSDY*DSDC
                H_[2,1] = H_[1,2]
                H_[2,2] = SPS*D2SDC2+2.0*DSDC*DSDC
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eISQ(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((1,2))
        IV_ = np.zeros(1)
        U_[0,0] = U_[0,0]+1
        U_[0,1] = U_[0,1]-1
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

