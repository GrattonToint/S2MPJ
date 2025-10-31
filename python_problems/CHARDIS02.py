from s2mpjlib import *
class  CHARDIS02(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CHARDIS02
#    *********
# 
#    Distribution of (equal)charges on [-R,R]x[-R,R] (2D)
# 
#    SIF input: R. Felkel, Jun 1999.
#               correction by S. Gratton & Ph. Toint, May 2024
#    modifield version of CHARDIS0 (formulation corrected)
# 
#    classification = "C-COBR2-AY-V-V"
# 
#    Number of positive (or negative) charges -> Number of variables 2*NP1
# 
#           Alternative values for the SIF file parameters:
# IE NP1                 5              $-PARAMETER
# IE NP1                 9              $-PARAMETER
# IE NP1                 20             $-PARAMETER
# IE NP1                 30             $-PARAMETER
# IE NP1                 50             $-PARAMETER     original value
# IE NP1                 100            $-PARAMETER
# IE NP1                 200            $-PARAMETER
# IE NP1                 500            $-PARAMETER
# IE NP1                 1000           $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'CHARDIS02'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['NP1'] = int(20);  #  SIF file default value
        else:
            v_['NP1'] = int(args[0])
# IE NP1                 2000           $-PARAMETER
# IE NP1                 5000           $-PARAMETER
        v_['R'] = 10.0
        v_['R-'] = -10.0
        v_['N'] = -1+v_['NP1']
        v_['NReal'] = float(v_['N'])
        v_['NP1Real'] = float(v_['NP1'])
        v_['halfPI'] = np.arcsin(1.0)
        v_['PI'] = 2.0*v_['halfPI']
        v_['2PI'] = 4.0*v_['halfPI']
        v_['4PI'] = 8.0*v_['halfPI']
        v_['4PIqN'] = v_['4PI']/v_['NReal']
        v_['2PIqN'] = v_['2PI']/v_['NReal']
        v_['PIqN'] = v_['PI']/v_['NReal']
        v_['RqN'] = v_['R']/v_['NReal']
        v_['1'] = 1
        v_['2'] = 2
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['NP1'])+1):
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
        for I in range(int(v_['1']),int(v_['NP1'])+1):
            v_['I+'] = 1+I
            for J in range(int(v_['I+']),int(v_['NP1'])+1):
                [ig,ig_,_] = s2mpj_ii('O'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
                self.gscale = arrset(self.gscale,ig,float(0.01))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        for I in range(int(v_['1']),int(v_['NP1'])+1):
            self.xlower[ix_['X'+str(I)]] = v_['R-']
            self.xupper[ix_['X'+str(I)]] = v_['R']
            self.xlower[ix_['Y'+str(I)]] = v_['R-']
            self.xupper[ix_['Y'+str(I)]] = v_['R']
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        for I in range(int(v_['1']),int(v_['NP1'])+1):
            v_['RealI-'] = float(I)
            v_['RealNP1-I'] = v_['NP1Real']-v_['RealI-']
            v_['PHII-'] = v_['2PIqN']*v_['RealI-']
            v_['RI-'] = v_['RqN']*v_['RealNP1-I']
            v_['XSTT'] = np.cos(v_['PHII-'])
            v_['YSTT'] = np.sin(v_['PHII-'])
            v_['XST'] = v_['XSTT']*v_['RI-']
            v_['YST'] = v_['YSTT']*v_['RI-']
            v_['XS'] = 0.5*v_['XST']
            v_['YS'] = 0.5*v_['YST']
            self.x0[ix_['X'+str(I)]] = float(v_['XS'])
            self.x0[ix_['Y'+str(I)]] = float(v_['YS'])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eDIFSQR', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['1']),int(v_['NP1'])+1):
            v_['I+'] = 1+I
            for J in range(int(v_['I+']),int(v_['NP1'])+1):
                ename = 'X'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eDIFSQR')
                ielftype = arrset(ielftype,ie,iet_["eDIFSQR"])
                vname = 'X'+str(I)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'X'+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                ename = 'Y'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eDIFSQR')
                ielftype = arrset(ielftype,ie,iet_["eDIFSQR"])
                vname = 'Y'+str(I)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'Y'+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gREZIP',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['NP1'])+1):
            v_['I+'] = 1+I
            for J in range(int(v_['I+']),int(v_['NP1'])+1):
                ig = ig_['O'+str(I)+','+str(J)]
                self.grftype = arrset(self.grftype,ig,'gREZIP')
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['X'+str(I)+','+str(J)])
                self.grelw = loaset(self.grelw,ig,posel,1.)
                posel = posel+1
                self.grelt = loaset(self.grelt,ig,posel,ie_['Y'+str(I)+','+str(J)])
                self.grelw = loaset(self.grelw,ig,posel, 1.)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass   = "C-COBR2-AY-V-V"
        self.objderlvl = 2


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eDIFSQR(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = (EV_[0]-EV_[1])*(EV_[0]-EV_[1])
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*(EV_[0]-EV_[1])
            g_[1] = -2.0*(EV_[0]-EV_[1])
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 2.0
                H_[0,1] = -2.0
                H_[1,0] = H_[0,1]
                H_[1,1] = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gREZIP(self,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= 1.0/GVAR_
        if nargout>1:
            g_ = -1.0/(GVAR_*GVAR_)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 2.0/(GVAR_*GVAR_*GVAR_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

