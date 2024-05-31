from s2xlib import *
class  OPTCNTRL(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : OPTCNTRL
#    *********
#    An optimal control problem
# 
#    Source:
#    B. Murtagh and M. Saunders,
#    Mathematical Programming studies 16, pp 84-117,
#    (example 5.11)
# 
#    SIF input: Nick Gould, June 1990.
# 
#    classification = "QQR2-AN-32-20"
# 
#    useful parameters
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'OPTCNTRL'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'OPTCNTRL'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['T'] = 10
        v_['T-1'] = -1+v_['T']
        v_['0'] = 0
        v_['1'] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for t in range(int(v_['0']),int(v_['T'])+1):
            [iv,ix_,_] = s2x_ii('x'+str(t),ix_)
            pb.xnames=arrset(pb.xnames,iv,'x'+str(t))
            [iv,ix_,_] = s2x_ii('y'+str(t),ix_)
            pb.xnames=arrset(pb.xnames,iv,'y'+str(t))
        for t in range(int(v_['0']),int(v_['T-1'])+1):
            [iv,ix_,_] = s2x_ii('u'+str(t),ix_)
            pb.xnames=arrset(pb.xnames,iv,'u'+str(t))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2x_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        for t in range(int(v_['0']),int(v_['T-1'])+1):
            v_['t+1'] = 1+t
            [ig,ig_,_] = s2x_ii('B'+str(t),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'B'+str(t))
            iv = ix_['x'+str(int(v_['t+1']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['x'+str(t)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            iv = ix_['y'+str(t)]
            pbm.A[ig,iv] = float(-0.2)+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('C'+str(t),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'C'+str(t))
            iv = ix_['y'+str(int(v_['t+1']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['y'+str(t)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            iv = ix_['x'+str(t)]
            pbm.A[ig,iv] = float(0.004)+pbm.A[ig,iv]
            iv = ix_['u'+str(t)]
            pbm.A[ig,iv] = float(-0.2)+pbm.A[ig,iv]
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
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('inf'))
        pb.xupper = np.full((pb.n,1),float('inf'))
        for t in range(int(v_['0']),int(v_['T-1'])+1):
            pb.xlower[ix_['x'+str(t)]] = -float('Inf')
            pb.xupper[ix_['x'+str(t)]] = +float('Inf')
            pb.xlower[ix_['y'+str(t)]] = -1.0
            pb.xlower[ix_['u'+str(t)]] = -0.2
            pb.xupper[ix_['u'+str(t)]] = 0.2
        pb.xlower[ix_['x'+str(int(v_['0']))]] = 10.0
        pb.xupper[ix_['x'+str(int(v_['0']))]] = 10.0
        pb.xlower[ix_['y'+str(int(v_['0']))]] = 0.0
        pb.xupper[ix_['y'+str(int(v_['0']))]] = 0.0
        pb.xlower[ix_['y'+str(int(v_['T']))]] = 0.0
        pb.xupper[ix_['y'+str(int(v_['T']))]] = 0.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        for t in range(int(v_['1']),int(v_['T-1'])+1):
            if('y'+str(t) in ix_):
                pb.x0[ix_['y'+str(t)]] = float(-1.0)
            else:
                pb.y0  = (
                      arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['y'+str(t)]),float(-1.0)))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eSQR', iet_)
        elftv = loaset(elftv,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for t in range(int(v_['0']),int(v_['T'])+1):
            ename = 'o'+str(t)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eSQR')
            ielftype = arrset(ielftype, ie, iet_["eSQR"])
            vname = 'x'+str(t)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for t in range(int(v_['0']),int(v_['T-1'])+1):
            ename = 'c'+str(t)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eSQR')
            ielftype = arrset(ielftype, ie, iet_["eSQR"])
            vname = 'y'+str(t)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for t in range(int(v_['0']),int(v_['T'])+1):
            ig = ig_['OBJ']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['o'+str(t)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.5))
        for t in range(int(v_['0']),int(v_['T-1'])+1):
            ig = ig_['C'+str(t)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['c'+str(t)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.01))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "QQR2-AN-32-20"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSQR(pbm,nargout,*args):

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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

