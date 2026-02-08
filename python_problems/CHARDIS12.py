from s2mpjlib import *
class  CHARDIS12(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CHARDIS12
#    *********
# 
#    Distribution of charges on a round plate (2D)
# 
#    SIF input: R. Felkel, Jun 1999.
#               correction by S. Gratton & Ph. Toint, May 2024
#    modifield version of CHARDIS1 (formulation corrected)
# 
#    classification = "C-COQR2-AY-V-V"
# 
#    Number of positive (or negative) charges -> Number of variables 2*NP1
# 
#           Alternative values for the SIF file parameters:
# IE NP1                 5              $-PARAMETER
# IE NP1                 8              $-PARAMETER
# IE NP1                 20             $-PARAMETER     original value
# IE NP1                 50             $-PARAMETER
# IE NP1                 100            $-PARAMETER
# IE NP1                 200            $-PARAMETER
# IE NP1                 500            $-PARAMETER
# IE NP1                 1000           $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'CHARDIS12'

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
        v_['R'] = 1.0
        v_['R2'] = v_['R']*v_['R']
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
        for I in range(int(v_['2']),int(v_['NP1'])+1):
            [ig,ig_,_] = s2mpj_ii('RES'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'RES'+str(I))
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
        for I in range(int(v_['2']),int(v_['NP1'])+1):
            self.gconst = arrset(self.gconst,ig_['RES'+str(I)],float(v_['R2']))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        for I in range(int(v_['2']),int(v_['NP1'])+1):
            self.xlower[ix_['X'+str(I)]] = -float('Inf')
            self.xupper[ix_['X'+str(I)]] = +float('Inf')
            self.xlower[ix_['Y'+str(I)]] = -float('Inf')
            self.xupper[ix_['Y'+str(I)]] = +float('Inf')
        self.xlower[ix_['X1']] = v_['R']
        self.xupper[ix_['X1']] = v_['R']
        self.xlower[ix_['Y1']] = 0.0
        self.xupper[ix_['Y1']] = 0.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        for I in range(int(v_['2']),int(v_['NP1'])+1):
            v_['I-'] = -1+I
            v_['RealI-'] = float(v_['I-'])
            v_['RealNP1-I'] = v_['NP1Real']-v_['RealI-']
            v_['PHII-'] = v_['2PIqN']*v_['RealI-']
            v_['RI-'] = v_['RqN']*v_['RealNP1-I']
            v_['XST'] = np.cos(v_['PHII-'])
            v_['YST'] = np.sin(v_['PHII-'])
            v_['XS'] = v_['XST']*v_['RI-']
            v_['YS'] = v_['YST']*v_['RI-']
            self.x0[ix_['X'+str(I)]] = float(v_['XS'])
            self.x0[ix_['Y'+str(I)]] = float(v_['YS'])
        self.x0[ix_['X1']] = float(v_['R'])
        self.x0[ix_['Y1']] = float(0.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQR', iet_)
        elftv = loaset(elftv,it,0,'V1')
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
        for I in range(int(v_['2']),int(v_['NP1'])+1):
            ename = 'RX'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_["eSQR"])
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'RY'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_["eSQR"])
            vname = 'Y'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
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
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,1.)
                posel = posel+1
                self.grelt = loaset(self.grelt,ig,posel,ie_['Y'+str(I)+','+str(J)])
                self.grelw = loaset(self.grelw,ig,posel, 1.)
        for I in range(int(v_['2']),int(v_['NP1'])+1):
            ig = ig_['RES'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['RX'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(1.0))
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['RY'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.cupper[np.arange(self.nle)] = np.zeros((self.nle,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-COQR2-AY-V-V"
        self.objderlvl = 2
        self.conderlvl = [2]


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSQR(self, nargout,*args):

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
    def eDIFSQR(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = (EV_[0,0]-EV_[1,0])*(EV_[0,0]-EV_[1,0])
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*(EV_[0,0]-EV_[1,0])
            g_[1] = -2.0*(EV_[0,0]-EV_[1,0])
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

