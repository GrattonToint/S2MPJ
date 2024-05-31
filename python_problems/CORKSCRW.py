from s2xlib import *
class  CORKSCRW(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CORKSCRW
#    *********
# 
#    A nonlinear optimal control problem with both state- and
#    control constraints.
#    The problem is to control (using an applied force of limited
#    magnitude) a mass moving in the 3D space, such that its
#    trajectory lies within a prescribed distance TOL of the
#    corkscreww-like curve defined by
#               y = sin(x), z = cos(x),
#    and such that it stops at a given abscissa XT in minimum time.
#    The mass is initially stationary at (0,0,1).
# 
#    Source:
#    Ph. Toint, private communication.
# 
#    SIF input: Ph. Toint, April 1991.
# 
#    classification = "SOR2-AN-V-V"
# 
#    Number of time intervals
#    The number of variables is 9T+6, of which 9 are fixed.
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'CORKSCRW'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'CORKSCRW'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['T'] = int(10);  #  SIF file default value
        else:
            v_['T'] = int(args[0])
#           Alternative values for the SIF file parameters:
# IE T                   50             $-PARAMETER n = 456
# IE T                   100            $-PARAMETER n = 906
# IE T                   500            $-PARAMETER n = 4506
# IE T                   1000           $-PARAMETER n = 9006
        if nargin<2:
            v_['XT'] = float(10.0);  #  SIF file default value
        else:
            v_['XT'] = float(args[1])
        if nargin<3:
            v_['MASS'] = float(0.37);  #  SIF file default value
        else:
            v_['MASS'] = float(args[2])
        if nargin<4:
            v_['TOL'] = float(0.1);  #  SIF file default value
        else:
            v_['TOL'] = float(args[3])
        v_['0'] = 0
        v_['1'] = 1
        v_['RT'] = float(v_['T'])
        v_['T+1'] = 1.0+v_['RT']
        v_['H'] = v_['XT']/v_['RT']
        v_['1/H'] = 1.0/v_['H']
        v_['-1/H'] = -1.0*v_['1/H']
        v_['M/H'] = v_['MASS']/v_['H']
        v_['-M/H'] = -1.0*v_['M/H']
        v_['TOLSQ'] = v_['TOL']*v_['TOL']
        v_['XTT+1'] = v_['XT']*v_['T+1']
        v_['W'] = 0.5*v_['XTT+1']
        for I in range(int(v_['1']),int(v_['T'])+1):
            v_['RI'] = float(I)
            v_['TI'] = v_['RI']*v_['H']
            v_['W/T'+str(I)] = v_['W']/v_['TI']
        v_['FMAX'] = v_['XT']/v_['RT']
        v_['-FMAX'] = -1.0*v_['FMAX']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['0']),int(v_['T'])+1):
            [iv,ix_,_] = s2x_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
            [iv,ix_,_] = s2x_ii('Y'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'Y'+str(I))
            [iv,ix_,_] = s2x_ii('Z'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'Z'+str(I))
            [iv,ix_,_] = s2x_ii('VX'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'VX'+str(I))
            [iv,ix_,_] = s2x_ii('VY'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'VY'+str(I))
            [iv,ix_,_] = s2x_ii('VZ'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'VZ'+str(I))
        for I in range(int(v_['1']),int(v_['T'])+1):
            [iv,ix_,_] = s2x_ii('UX'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'UX'+str(I))
            [iv,ix_,_] = s2x_ii('UY'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'UY'+str(I))
            [iv,ix_,_] = s2x_ii('UZ'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'UZ'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['T'])+1):
            [ig,ig_,_] = s2x_ii('OX'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            pbm.gscale = arrset(pbm.gscale,ig,float(v_['W/T'+str(I)]))
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['T'])+1):
            v_['I-1'] = -1+I
            [ig,ig_,_] = s2x_ii('ACX'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'ACX'+str(I))
            iv = ix_['VX'+str(I)]
            pbm.A[ig,iv] = float(v_['M/H'])+pbm.A[ig,iv]
            iv = ix_['VX'+str(int(v_['I-1']))]
            pbm.A[ig,iv] = float(v_['-M/H'])+pbm.A[ig,iv]
            iv = ix_['UX'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('ACY'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'ACY'+str(I))
            iv = ix_['VY'+str(I)]
            pbm.A[ig,iv] = float(v_['M/H'])+pbm.A[ig,iv]
            iv = ix_['VY'+str(int(v_['I-1']))]
            pbm.A[ig,iv] = float(v_['-M/H'])+pbm.A[ig,iv]
            iv = ix_['UY'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('ACZ'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'ACZ'+str(I))
            iv = ix_['VZ'+str(I)]
            pbm.A[ig,iv] = float(v_['M/H'])+pbm.A[ig,iv]
            iv = ix_['VZ'+str(int(v_['I-1']))]
            pbm.A[ig,iv] = float(v_['-M/H'])+pbm.A[ig,iv]
            iv = ix_['UZ'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('PSX'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'PSX'+str(I))
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = float(v_['1/H'])+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['I-1']))]
            pbm.A[ig,iv] = float(v_['-1/H'])+pbm.A[ig,iv]
            iv = ix_['VX'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('PSY'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'PSY'+str(I))
            iv = ix_['Y'+str(I)]
            pbm.A[ig,iv] = float(v_['1/H'])+pbm.A[ig,iv]
            iv = ix_['Y'+str(int(v_['I-1']))]
            pbm.A[ig,iv] = float(v_['-1/H'])+pbm.A[ig,iv]
            iv = ix_['VY'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('PSZ'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'PSZ'+str(I))
            iv = ix_['Z'+str(I)]
            pbm.A[ig,iv] = float(v_['1/H'])+pbm.A[ig,iv]
            iv = ix_['Z'+str(int(v_['I-1']))]
            pbm.A[ig,iv] = float(v_['-1/H'])+pbm.A[ig,iv]
            iv = ix_['VZ'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('SC'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'SC'+str(I))
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
        for I in range(int(v_['1']),int(v_['T'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['OX'+str(I)],float(v_['XT']))
            pbm.gconst = arrset(pbm.gconst,ig_['SC'+str(I)],float(v_['TOLSQ']))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower[ix_['X'+str(int(v_['0']))]] = 0.0
        pb.xupper[ix_['X'+str(int(v_['0']))]] = 0.0
        pb.xlower[ix_['Y'+str(int(v_['0']))]] = 0.0
        pb.xupper[ix_['Y'+str(int(v_['0']))]] = 0.0
        pb.xlower[ix_['Z'+str(int(v_['0']))]] = 1.0
        pb.xupper[ix_['Z'+str(int(v_['0']))]] = 1.0
        pb.xlower[ix_['VX'+str(int(v_['0']))]] = 0.0
        pb.xupper[ix_['VX'+str(int(v_['0']))]] = 0.0
        pb.xlower[ix_['VY'+str(int(v_['0']))]] = 0.0
        pb.xupper[ix_['VY'+str(int(v_['0']))]] = 0.0
        pb.xlower[ix_['VZ'+str(int(v_['0']))]] = 0.0
        pb.xupper[ix_['VZ'+str(int(v_['0']))]] = 0.0
        pb.xlower[ix_['VX'+str(int(v_['T']))]] = 0.0
        pb.xupper[ix_['VX'+str(int(v_['T']))]] = 0.0
        pb.xlower[ix_['VY'+str(int(v_['T']))]] = 0.0
        pb.xupper[ix_['VY'+str(int(v_['T']))]] = 0.0
        pb.xlower[ix_['VZ'+str(int(v_['T']))]] = 0.0
        pb.xupper[ix_['VZ'+str(int(v_['T']))]] = 0.0
        for I in range(int(v_['1']),int(v_['T'])+1):
            pb.xlower[ix_['UX'+str(I)]] = v_['-FMAX']
            pb.xupper[ix_['UX'+str(I)]] = v_['FMAX']
            pb.xlower[ix_['UY'+str(I)]] = v_['-FMAX']
            pb.xupper[ix_['UY'+str(I)]] = v_['FMAX']
            pb.xlower[ix_['UZ'+str(I)]] = v_['-FMAX']
            pb.xupper[ix_['UZ'+str(I)]] = v_['FMAX']
            pb.xlower[ix_['X'+str(I)]] = 0.0
            pb.xupper[ix_['X'+str(I)]] = v_['XT']
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        pb.x0[ix_['X'+str(int(v_['0']))]] = float(0.0)
        pb.x0[ix_['Y'+str(int(v_['0']))]] = float(0.0)
        pb.x0[ix_['Z'+str(int(v_['0']))]] = float(1.0)
        pb.x0[ix_['VX'+str(int(v_['0']))]] = float(0.0)
        pb.x0[ix_['VY'+str(int(v_['0']))]] = float(0.0)
        pb.x0[ix_['VZ'+str(int(v_['0']))]] = float(0.0)
        pb.x0[ix_['VX'+str(int(v_['T']))]] = float(0.0)
        pb.x0[ix_['VY'+str(int(v_['T']))]] = float(0.0)
        pb.x0[ix_['VZ'+str(int(v_['T']))]] = float(0.0)
        for I in range(int(v_['1']),int(v_['T'])+1):
            v_['RI'] = float(I)
            v_['TI'] = v_['RI']*v_['H']
            if('X'+str(I) in ix_):
                pb.x0[ix_['X'+str(I)]] = float(v_['TI'])
            else:
                pb.y0  = (
                      arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X'+str(I)]),float(v_['TI'])))
            if('VX'+str(I) in ix_):
                pb.x0[ix_['VX'+str(I)]] = float(1.0)
            else:
                pb.y0  = (
                      arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['VX'+str(I)]),float(1.0)))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eERRSIN', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        [it,iet_,_] = s2x_ii( 'eERRCOS', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Z')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['1']),int(v_['T'])+1):
            ename = 'ES'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eERRSIN')
            ielftype = arrset(ielftype, ie, iet_["eERRSIN"])
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'Y'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'EC'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eERRCOS')
            ielftype = arrset(ielftype, ie, iet_["eERRCOS"])
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'Z'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2x_ii('gL2',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['T'])+1):
            ig = ig_['OX'+str(I)]
            pbm.grftype = arrset(pbm.grftype,ig,'gL2')
            ig = ig_['SC'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['ES'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EC'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.cupper[np.arange(pb.nle)] = np.zeros((pb.nle,1))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "SOR2-AN-V-V"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eERRSIN(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        SINX = np.sin(EV_[0])
        COSX = np.cos(EV_[0])
        ERR = EV_[1]-SINX
        f_   = ERR*ERR
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -2.0*ERR*COSX
            g_[1] = 2.0*ERR
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 2.0*(COSX**2+ERR*SINX)
                H_[0,1] = -2.0*COSX
                H_[1,0] = H_[0,1]
                H_[1,1] = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eERRCOS(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        SINX = np.sin(EV_[0])
        COSX = np.cos(EV_[0])
        ERR = EV_[1]-COSX
        f_   = ERR*ERR
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*ERR*SINX
            g_[1] = 2.0*ERR
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 2.0*(SINX**2+ERR*COSX)
                H_[0,1] = 2.0*SINX
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

