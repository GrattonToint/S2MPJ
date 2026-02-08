from s2mpjlib import *
class  AIRCRFTB(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : AIRCRFTB
#    *********
# 
#    The aircraft stability problem by Rheinboldt, as a function
#    of the elevator, aileron and rudder deflection controls.
# 
#    Source: problem 9 in
#    J.J. More',"A collection of nonlinear model problems"
#    Proceedings of the AMS-SIAM Summer Seminar on the Computational
#    Solution of Nonlinear Systems of Equations, Colorado, 1988.
#    Argonne National Laboratory MCS-P60-0289, 1989.
# 
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "C-CSXR2-RN-8-0"
# 
#    Values for the controls
#    1) Elevator
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'AIRCRFTB'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['ELVVAL'] = -0.05
        v_['AILVAL'] = 0.1
        v_['RUDVAL'] = 0.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [iv,ix_,_] = s2mpj_ii('ROLLRATE',ix_)
        self.xnames=arrset(self.xnames,iv,'ROLLRATE')
        [iv,ix_,_] = s2mpj_ii('PITCHRAT',ix_)
        self.xnames=arrset(self.xnames,iv,'PITCHRAT')
        [iv,ix_,_] = s2mpj_ii('YAWRATE',ix_)
        self.xnames=arrset(self.xnames,iv,'YAWRATE')
        [iv,ix_,_] = s2mpj_ii('ATTCKANG',ix_)
        self.xnames=arrset(self.xnames,iv,'ATTCKANG')
        [iv,ix_,_] = s2mpj_ii('SSLIPANG',ix_)
        self.xnames=arrset(self.xnames,iv,'SSLIPANG')
        [iv,ix_,_] = s2mpj_ii('ELEVATOR',ix_)
        self.xnames=arrset(self.xnames,iv,'ELEVATOR')
        [iv,ix_,_] = s2mpj_ii('AILERON',ix_)
        self.xnames=arrset(self.xnames,iv,'AILERON')
        [iv,ix_,_] = s2mpj_ii('RUDDERDF',ix_)
        self.xnames=arrset(self.xnames,iv,'RUDDERDF')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('G1',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['ROLLRATE']])
        valA = np.append(valA,float(-3.933))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['PITCHRAT']])
        valA = np.append(valA,float(0.107))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['YAWRATE']])
        valA = np.append(valA,float(0.126))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['SSLIPANG']])
        valA = np.append(valA,float(-9.99))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['AILERON']])
        valA = np.append(valA,float(-45.83))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['RUDDERDF']])
        valA = np.append(valA,float(-7.64))
        [ig,ig_,_] = s2mpj_ii('G2',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['PITCHRAT']])
        valA = np.append(valA,float(-0.987))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['ATTCKANG']])
        valA = np.append(valA,float(-22.95))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['ELEVATOR']])
        valA = np.append(valA,float(-28.37))
        [ig,ig_,_] = s2mpj_ii('G3',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['ROLLRATE']])
        valA = np.append(valA,float(0.002))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['YAWRATE']])
        valA = np.append(valA,float(-0.235))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['SSLIPANG']])
        valA = np.append(valA,float(5.67))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['AILERON']])
        valA = np.append(valA,float(-0.921))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['RUDDERDF']])
        valA = np.append(valA,float(-6.51))
        [ig,ig_,_] = s2mpj_ii('G4',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['PITCHRAT']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['ATTCKANG']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['ELEVATOR']])
        valA = np.append(valA,float(-1.168))
        [ig,ig_,_] = s2mpj_ii('G5',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['YAWRATE']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['SSLIPANG']])
        valA = np.append(valA,float(-0.196))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['AILERON']])
        valA = np.append(valA,float(-0.0071))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        selfnob      = ngrp
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        self.xlower[ix_['ELEVATOR']] = v_['ELVVAL']
        self.xupper[ix_['ELEVATOR']] = v_['ELVVAL']
        self.xlower[ix_['AILERON']] = v_['AILVAL']
        self.xupper[ix_['AILERON']] = v_['AILVAL']
        self.xlower[ix_['RUDDERDF']] = v_['RUDVAL']
        self.xupper[ix_['RUDDERDF']] = v_['RUDVAL']
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(0.0))
        self.x0[ix_['ELEVATOR']] = float(v_['ELVVAL'])
        self.x0[ix_['AILERON']] = float(v_['AILVAL'])
        self.x0[ix_['RUDDERDF']] = float(v_['RUDVAL'])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        ename = 'E1A'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'PITCHRAT'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'YAWRATE'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E1B'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'YAWRATE'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'ATTCKANG'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E1C'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'ATTCKANG'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'SSLIPANG'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E1D'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'PITCHRAT'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'ATTCKANG'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E2A'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'ROLLRATE'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'YAWRATE'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E2B'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'ROLLRATE'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'SSLIPANG'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E3A'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'ROLLRATE'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'PITCHRAT'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E3B'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'ROLLRATE'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'ATTCKANG'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E3C'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'PITCHRAT'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'ATTCKANG'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E4'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'ROLLRATE'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'SSLIPANG'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E5'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_['en2PR'])
        vname = 'ROLLRATE'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'ATTCKANG'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
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
        for ig in range(0,ngrp):
            self.grftype = arrset(self.grftype,ig,'gL2')
        ig = ig_['G1']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E1A'])
        self.grelw = loaset(self.grelw,ig,posel,float(-0.727))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E1B'])
        self.grelw = loaset(self.grelw,ig,posel,float(8.39))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E1C'])
        self.grelw = loaset(self.grelw,ig,posel,float(-684.4))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E1D'])
        self.grelw = loaset(self.grelw,ig,posel,float(63.5))
        ig = ig_['G2']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E2A'])
        self.grelw = loaset(self.grelw,ig,posel,float(0.949))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E2B'])
        self.grelw = loaset(self.grelw,ig,posel,float(0.173))
        ig = ig_['G3']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E3A'])
        self.grelw = loaset(self.grelw,ig,posel,float(-0.716))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E3B'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.578))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E3C'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.132))
        ig = ig_['G4']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E4'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        ig = ig_['G5']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E5'])
        self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 0.0
#    Solution
# LO SOLTN               6.4099D-02
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass   = "C-CSXR2-RN-8-0"
        self.objderlvl = 2

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

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

