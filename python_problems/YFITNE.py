from s2xlib import *
class  YFITNE(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    This problem arises in measuring angles and distances to a vibrating beam
#    using a laser-Doppler velocimeter.
#    This is a nonlinear equation variant of the bounded constrained
#    problem YFIT.
# 
#    Source:
#    an exercize for L. Watson course on LANCELOT in the Spring 1993.
# 
#    SIF input: B. E. Lindholm, Virginia Tech., Spring 1993,
#               modified by Ph. Toint, March 1994.
#               derivatives corrected by Nick Gould, June 2019.
# 
#    classification = "NOR2-MN-3-16"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'YFITNE'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'YFITNE'
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
        [iv,ix_,_] = s2x_ii('alpha',ix_)
        pb.xnames=arrset(pb.xnames,iv,'alpha')
        [iv,ix_,_] = s2x_ii('beta',ix_)
        pb.xnames=arrset(pb.xnames,iv,'beta')
        [iv,ix_,_] = s2x_ii('dist',ix_)
        pb.xnames=arrset(pb.xnames,iv,'dist')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for i in range(int(v_['zero']),int(v_['p'])+1):
            [ig,ig_,_] = s2x_ii('diff'+str(i),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'diff'+str(i))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        legrps = find(gtype,lambda x:x=='<=')
        eqgrps = find(gtype,lambda x:x=='==')
        gegrps = find(gtype,lambda x:x=='>=')
        pb.nle = len(legrps)
        pb.neq = len(eqgrps)
        pb.nge = len(gegrps)
        pb.m   = pb.nle+pb.neq+pb.nge
        pbm.congrps = find(gtype,lambda x:(x=='<=' or x=='==' or x=='>='))
        pb.cnames= cnames[pbm.congrps]
        pb.nob = ngrp-pb.m
        pbm.objgrps = find(gtype,lambda x:x=='<>')
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        for i in range(int(v_['zero']),int(v_['p'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['diff'+str(i)],float(v_['y'+str(i)]))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        pb.x0[ix_['alpha']] = float(0.60)
        pb.x0[ix_['beta']] = float(-0.60)
        pb.x0[ix_['dist']] = float(20.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'etanab', iet_)
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
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'etanab')
            ielftype = arrset(ielftype, ie, iet_["etanab"])
            vname = 'alpha'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='a1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'beta'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='b1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'dist'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='d1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='point')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['index']))
            posep = find(elftp[ielftype[ie]],lambda x:x=='count')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['realp']))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for i in range(int(v_['zero']),int(v_['p'])+1):
            ig = ig_['diff'+str(i)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['est'+str(i)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "NOR2-MN-3-16"
        self.pb = pb; self.pbm = pbm

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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

