from s2mpjlib import *
class  YFIT(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    A nonlinear least-squares problem.
# 
#    Source:
#    an exercize for L. Watson course on LANCELOT in the Spring 1993.
# 
#    SIF input: Brian E. Lindholm, Virginia Tech., Spring 1993.
#               derivatives corrected by Nick Gould, June 2019.
# 
#    classification = "SBR2-MN-3-0"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'YFIT'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['zero'] = 0
        v_['p'] = 16
        v_['realp'] = 16.0
        v_['y0'] = 21.158931
        v_['y1'] = 17.591719
        v_['y2'] = 14.046854
        v_['y3'] = 10.519732
        v_['y4'] = 7.0058392
        v_['y5'] = 3.5007293
        v_['y6'] = 0.0000000
        v_['y7'] = -3.5007293
        v_['y8'] = -7.0058392
        v_['y9'] = -10.519732
        v_['y10'] = -14.046854
        v_['y11'] = -17.591719
        v_['y12'] = -21.158931
        v_['y13'] = -24.753206
        v_['y14'] = -28.379405
        v_['y15'] = -32.042552
        v_['y16'] = -35.747869
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2mpj_ii('alpha',ix_)
        pb.xnames=arrset(pb.xnames,iv,'alpha')
        [iv,ix_,_] = s2mpj_ii('beta',ix_)
        pb.xnames=arrset(pb.xnames,iv,'beta')
        [iv,ix_,_] = s2mpj_ii('dist',ix_)
        pb.xnames=arrset(pb.xnames,iv,'dist')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for i in range(int(v_['zero']),int(v_['p'])+1):
            [ig,ig_,_] = s2mpj_ii('diff'+str(i),ig_)
            gtype = arrset(gtype,ig,'<>')
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        for i in range(int(v_['zero']),int(v_['p'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['diff'+str(i)],float(v_['y'+str(i)]))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),float('inf'))
        pb.xlower[ix_['alpha']] = -float('Inf')
        pb.xupper[ix_['alpha']] = +float('Inf')
        pb.xlower[ix_['beta']] = -float('Inf')
        pb.xupper[ix_['beta']] = +float('Inf')
        pb.xlower[ix_['dist']] = 0.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.x0[ix_['alpha']] = float(0.60)
        pb.x0[ix_['beta']] = float(-0.60)
        pb.x0[ix_['dist']] = float(20.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'etanab', iet_)
        elftv = loaset(elftv,it,0,'a1')
        elftv = loaset(elftv,it,1,'b1')
        elftv = loaset(elftv,it,2,'d1')
        elftp = []
        elftp = loaset(elftp,it,0,'point')
        elftp = loaset(elftp,it,1,'count')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for i in range(int(v_['zero']),int(v_['p'])+1):
            v_['index'] = float(i)
            ename = 'est'+str(i)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'etanab')
            ielftype = arrset(ielftype, ie, iet_["etanab"])
            vname = 'alpha'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='a1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'beta'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='b1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'dist'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='d1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='point')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['index']))
            posep = find(elftp[ielftype[ie]],lambda x:x=='count')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['realp']))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gsquare',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for i in range(int(v_['zero']),int(v_['p'])+1):
            ig = ig_['diff'+str(i)]
            pbm.grftype = arrset(pbm.grftype,ig,'gsquare')
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['est'+str(i)])
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
# LO SOLUTION            0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "SBR2-MN-3-0"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def etanab(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        frac = pbm.elpar[iel_][0]/pbm.elpar[iel_][1]
        ttan = np.tan(EV_[0]*(1.0-frac)+EV_[1]*frac)
        tsec = 1.0/np.cos(EV_[0]*(1.0-frac)+EV_[1]*frac)
        tsec2 = tsec*tsec
        f_   = EV_[2]*ttan
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[2]*(1.0-frac)*tsec2
            g_[1] = EV_[2]*frac*tsec2
            g_[2] = ttan
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = 2.0*EV_[2]*((1.0-frac)**2)*tsec2*ttan
                H_[1,1] = 2.0*EV_[2]*(frac**2)*tsec2*ttan
                H_[0,1] = 2.0*EV_[2]*(1.0-frac)*frac*tsec2*ttan
                H_[1,0] = H_[0,1]
                H_[0,2] = (1.0-frac)*tsec2
                H_[2,0] = H_[0,2]
                H_[1,2] = frac*tsec2
                H_[2,1] = H_[1,2]
                H_[2,2] = 0.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gsquare(pbm,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= GVAR_*GVAR_
        if nargout>1:
            g_ = 2.0*GVAR_
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

