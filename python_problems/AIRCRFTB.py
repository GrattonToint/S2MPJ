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
#    classification = "SXR2-RN-8-0"
# 
#    Values for the controls
#    1) Elevator
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'AIRCRFTB'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'AIRCRFTB'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['ELVVAL'] = -0.05
        v_['AILVAL'] = 0.1
        v_['RUDVAL'] = 0.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2mpj_ii('ROLLRATE',ix_)
        pb.xnames=arrset(pb.xnames,iv,'ROLLRATE')
        [iv,ix_,_] = s2mpj_ii('PITCHRAT',ix_)
        pb.xnames=arrset(pb.xnames,iv,'PITCHRAT')
        [iv,ix_,_] = s2mpj_ii('YAWRATE',ix_)
        pb.xnames=arrset(pb.xnames,iv,'YAWRATE')
        [iv,ix_,_] = s2mpj_ii('ATTCKANG',ix_)
        pb.xnames=arrset(pb.xnames,iv,'ATTCKANG')
        [iv,ix_,_] = s2mpj_ii('SSLIPANG',ix_)
        pb.xnames=arrset(pb.xnames,iv,'SSLIPANG')
        [iv,ix_,_] = s2mpj_ii('ELEVATOR',ix_)
        pb.xnames=arrset(pb.xnames,iv,'ELEVATOR')
        [iv,ix_,_] = s2mpj_ii('AILERON',ix_)
        pb.xnames=arrset(pb.xnames,iv,'AILERON')
        [iv,ix_,_] = s2mpj_ii('RUDDERDF',ix_)
        pb.xnames=arrset(pb.xnames,iv,'RUDDERDF')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('G1',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['ROLLRATE']
        pbm.A[ig,iv] = float(-3.933)+pbm.A[ig,iv]
        iv = ix_['PITCHRAT']
        pbm.A[ig,iv] = float(0.107)+pbm.A[ig,iv]
        iv = ix_['YAWRATE']
        pbm.A[ig,iv] = float(0.126)+pbm.A[ig,iv]
        iv = ix_['SSLIPANG']
        pbm.A[ig,iv] = float(-9.99)+pbm.A[ig,iv]
        iv = ix_['AILERON']
        pbm.A[ig,iv] = float(-45.83)+pbm.A[ig,iv]
        iv = ix_['RUDDERDF']
        pbm.A[ig,iv] = float(-7.64)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('G2',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['PITCHRAT']
        pbm.A[ig,iv] = float(-0.987)+pbm.A[ig,iv]
        iv = ix_['ATTCKANG']
        pbm.A[ig,iv] = float(-22.95)+pbm.A[ig,iv]
        iv = ix_['ELEVATOR']
        pbm.A[ig,iv] = float(-28.37)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('G3',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['ROLLRATE']
        pbm.A[ig,iv] = float(0.002)+pbm.A[ig,iv]
        iv = ix_['YAWRATE']
        pbm.A[ig,iv] = float(-0.235)+pbm.A[ig,iv]
        iv = ix_['SSLIPANG']
        pbm.A[ig,iv] = float(5.67)+pbm.A[ig,iv]
        iv = ix_['AILERON']
        pbm.A[ig,iv] = float(-0.921)+pbm.A[ig,iv]
        iv = ix_['RUDDERDF']
        pbm.A[ig,iv] = float(-6.51)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('G4',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['PITCHRAT']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['ATTCKANG']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['ELEVATOR']
        pbm.A[ig,iv] = float(-1.168)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('G5',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['YAWRATE']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['SSLIPANG']
        pbm.A[ig,iv] = float(-0.196)+pbm.A[ig,iv]
        iv = ix_['AILERON']
        pbm.A[ig,iv] = float(-0.0071)+pbm.A[ig,iv]
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower[ix_['ELEVATOR']] = v_['ELVVAL']
        pb.xupper[ix_['ELEVATOR']] = v_['ELVVAL']
        pb.xlower[ix_['AILERON']] = v_['AILVAL']
        pb.xupper[ix_['AILERON']] = v_['AILVAL']
        pb.xlower[ix_['RUDDERDF']] = v_['RUDVAL']
        pb.xupper[ix_['RUDDERDF']] = v_['RUDVAL']
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(0.0))
        pb.x0[ix_['ELEVATOR']] = float(v_['ELVVAL'])
        pb.x0[ix_['AILERON']] = float(v_['AILVAL'])
        pb.x0[ix_['RUDDERDF']] = float(v_['RUDVAL'])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        ename = 'E1A'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'PITCHRAT'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'YAWRATE'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E1B'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'YAWRATE'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'ATTCKANG'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E1C'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'ATTCKANG'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'SSLIPANG'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E1D'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'PITCHRAT'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'ATTCKANG'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E2A'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'ROLLRATE'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'YAWRATE'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E2B'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'ROLLRATE'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'SSLIPANG'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E3A'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'ROLLRATE'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'PITCHRAT'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E3B'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'ROLLRATE'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'ATTCKANG'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E3C'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'PITCHRAT'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'ATTCKANG'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E4'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'ROLLRATE'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'SSLIPANG'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E5'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'ROLLRATE'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'ATTCKANG'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gL2',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for ig in range(0,ngrp):
            pbm.grftype = arrset(pbm.grftype,ig,'gL2')
        ig = ig_['G1']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E1A'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-0.727))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E1B'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(8.39))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E1C'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-684.4))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E1D'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(63.5))
        ig = ig_['G2']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E2A'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.949))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E2B'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.173))
        ig = ig_['G3']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E3A'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-0.716))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E3B'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.578))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E3C'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.132))
        ig = ig_['G4']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E4'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        ig = ig_['G5']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E5'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN               6.4099D-02
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "SXR2-RN-8-0"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def en2PR(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]
            g_[1] = EV_[0]
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
    def gL2(pbm,nargout,*args):

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

