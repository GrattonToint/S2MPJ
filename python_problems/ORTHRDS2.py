from s2mpjlib import *
class  ORTHRDS2(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : ORTHRDS2
#    *********
# 
#    An orthogonal regression problem.
# 
#    The problem is to fit (orthogonally) a planar curve to a set of points
#    in the plane. This set of points is generated by perturbing a
#    first set lying exactly on the predefined curve.
#    The curve is referred to as a cardioid in the original paper,
#    but is in fact a circle.
# 
#    This problem is a modification of ORTHREGD.SIF.
#    The start coordinates of certain data points were moved to
#    force convergence to the center of the circle, where
#    the Jacobian of the constraints is singular.
# 
#    Source: adapted from:
#    M. Gulliksson,
#    "Algorithms for nonlinear Least-squares with Applications to
#    Orthogonal Regression",
#    UMINF-178.90, University of Umea, Sweden, 1990.
# 
#    SIF input: Ph. Toint, Mar 1991.
#               modified by T Plantagena, May 1994.
# 
#    classification = "C-CQOR2-AY-V-V"
# 
#    Number of data points
#    (number of variables = 2 NPTS + 3 )
# 
#           Alternative values for the SIF file parameters:
# IE NPTS                10             $-PARAMETER n = 23
# IE NPTS                50             $-PARAMETER n = 103
# IE NPTS                76             $-PARAMETER n = 155
# IE NPTS                100            $-PARAMETER n = 203    original value
# IE NPTS                250            $-PARAMETER n = 503
# IE NPTS                500            $-PARAMETER n = 1003
# IE NPTS                2500           $-PARAMETER n = 5003
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 25 XI 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'ORTHRDS2'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['NPTS'] = int(10);  #  SIF file default value
        else:
            v_['NPTS'] = int(args[0])
# IE NPTS                5000           $-PARAMETER n = 10003
        v_['TZ3'] = 1.7
        v_['PSEED'] = 237.1531
        v_['PSIZE'] = 0.2
        v_['1'] = 1
        v_['0'] = 0
        v_['TDPulo'] = 5
        v_['TDPuhi'] = 5
        v_['PI'] = 3.1415926535
        v_['2PI'] = 2.0*v_['PI']
        v_['RNPTS'] = float(v_['NPTS'])
        v_['ICR0'] = 1.0/v_['RNPTS']
        v_['INCR'] = v_['ICR0']*v_['2PI']
        v_['Z3SQ'] = v_['TZ3']*v_['TZ3']
        v_['1+TZ3SQ'] = 1.0+v_['Z3SQ']
        for I in range(int(v_['1']),int(v_['NPTS'])+1):
            v_['I-1'] = -1+I
            v_['RI-1'] = float(v_['I-1'])
            v_['THETA'] = v_['RI-1']*v_['INCR']
            v_['ST'] = np.sin(v_['THETA'])
            v_['CT'] = np.cos(v_['THETA'])
            v_['FACT'] = v_['1+TZ3SQ']+v_['CT']
            v_['R1'] = v_['FACT']*v_['CT']
            v_['R2'] = v_['FACT']*v_['ST']
            v_['XSEED'] = v_['THETA']*v_['PSEED']
            v_['SSEED'] = np.cos(v_['XSEED'])
            v_['PER-1'] = v_['PSIZE']*v_['SSEED']
            v_['PERT'] = 1.0+v_['PER-1']
            v_['XD'+str(I)] = v_['R1']*v_['PERT']
            v_['YD'+str(I)] = v_['R2']*v_['PERT']
        for I in range(int(v_['TDPulo']),int(v_['TDPuhi'])+1):
            v_['XD'+str(I)] = 1.1
            v_['YD'+str(I)] = 0.1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [iv,ix_,_] = s2mpj_ii('Z1',ix_)
        self.xnames=arrset(self.xnames,iv,'Z1')
        [iv,ix_,_] = s2mpj_ii('Z2',ix_)
        self.xnames=arrset(self.xnames,iv,'Z2')
        [iv,ix_,_] = s2mpj_ii('Z3',ix_)
        self.xnames=arrset(self.xnames,iv,'Z3')
        for I in range(int(v_['1']),int(v_['NPTS'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I))
            [iv,ix_,_] = s2mpj_ii('Y'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'Y'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['1']),int(v_['NPTS'])+1):
            [ig,ig_,_] = s2mpj_ii('OX'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(I)]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('OY'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['Y'+str(I)]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('E'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'E'+str(I))
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
            self.gconst = arrset(self.gconst,ig_['OX'+str(I)],float(v_['XD'+str(I)]))
            self.gconst = arrset(self.gconst,ig_['OY'+str(I)],float(v_['YD'+str(I)]))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        if('Z1' in ix_):
            self.x0[ix_['Z1']] = float(1.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Z1']),float(1.0)))
        if('Z2' in ix_):
            self.x0[ix_['Z2']] = float(0.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Z2']),float(0.0)))
        if('Z3' in ix_):
            self.x0[ix_['Z3']] = float(1.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Z3']),float(1.0)))
        for I in range(int(v_['1']),int(v_['NPTS'])+1):
            if('X'+str(I) in ix_):
                self.x0[ix_['X'+str(I)]] = float(v_['XD'+str(I)])
            else:
                self.y0  = (
                      arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X'+str(I)]),float(v_['XD'+str(I)])))
            if('Y'+str(I) in ix_):
                self.x0[ix_['Y'+str(I)]] = float(v_['YD'+str(I)])
            else:
                self.y0  = (
                      arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Y'+str(I)]),float(v_['YD'+str(I)])))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eTA', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        elftv = loaset(elftv,it,2,'ZA')
        elftv = loaset(elftv,it,3,'ZB')
        [it,iet_,_] = s2mpj_ii( 'eTB', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        elftv = loaset(elftv,it,2,'ZA')
        elftv = loaset(elftv,it,3,'ZB')
        elftv = loaset(elftv,it,4,'ZC')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['1']),int(v_['NPTS'])+1):
            ename = 'EA'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eTA')
            ielftype = arrset(ielftype,ie,iet_["eTA"])
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'Y'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='Y')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'Z1'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='ZA')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'Z2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='ZB')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'EB'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eTB')
            ielftype = arrset(ielftype,ie,iet_["eTB"])
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'Y'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='Y')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'Z1'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='ZA')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'Z2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='ZB')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'Z3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='ZC')[0]
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
            ig = ig_['OX'+str(I)]
            self.grftype = arrset(self.grftype,ig,'gL2')
            ig = ig_['OY'+str(I)]
            self.grftype = arrset(self.grftype,ig,'gL2')
            ig = ig_['E'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['EA'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = posel+1
            self.grelt = loaset(self.grelt,ig,posel,ie_['EB'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 0.0
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
        U_ = np.zeros((2,4))
        IV_ = np.zeros(2)
        U_[0,0] = U_[0,0]+1
        U_[0,2] = U_[0,2]-1
        U_[1,1] = U_[1,1]+1
        U_[1,3] = U_[1,3]-1
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        T = IV_[0]*IV_[0]+IV_[1]*IV_[1]
        f_   = T*T
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 4.0*T*IV_[0]
            g_[1] = 4.0*T*IV_[1]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 4.0*(T+2.0*IV_[0]*IV_[0])
                H_[0,1] = 8.0*IV_[0]*IV_[1]
                H_[1,0] = H_[0,1]
                H_[1,1] = 4.0*(T+2.0*IV_[1]*IV_[1])
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eTB(self, nargout,*args):

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
        T = IV_[0]*IV_[0]+IV_[1]*IV_[1]
        ZZSQ = IV_[2]*IV_[2]
        T1 = 1.0+ZZSQ
        T1SQ = T1*T1
        f_   = T*T1SQ
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*IV_[0]*T1SQ
            g_[1] = 2.0*IV_[1]*T1SQ
            g_[2] = 4.0*T*T1*IV_[2]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = 2.0*T1SQ
                H_[0,2] = 8.0*IV_[0]*T1*IV_[2]
                H_[2,0] = H_[0,2]
                H_[1,1] = 2.0*T1SQ
                H_[1,2] = 8.0*IV_[1]*T1*IV_[2]
                H_[2,1] = H_[1,2]
                H_[2,2] = 4.0*T*(2.0*ZZSQ+T1)
                H_ = U_.T.dot(H_).dot(U_)
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

