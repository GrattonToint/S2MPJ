from s2xlib import *
class  CLNLBEAM(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CLNLBEAM
#    *********
# 
#    An optimal control version of the CLamped NonLinear BEAM problem.
#    The energy of a beam of length 1 compressed by a force P is to be
#    minimized.  The control variable is the derivative of the deflection angle.
# 
#    The problem is discretized using the trapezoidal rule. It is non-convex.
# 
#    Source:
#    H. Maurer and H.D. Mittelman,
#    "The non-linear beam via optimal control with bound state variables",
#    Optimal Control Applications and Methods 12, pp. 19-31, 1991.
# 
#    SIF input: Ph. Toint, Nov 1993.
# 
#    classification = "OOR2-MN-V-V"
# 
#    Discretization: specify the number of interior points + 1
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'CLNLBEAM'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'CLNLBEAM'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['NI'] = int(10);  #  SIF file default value
        else:
            v_['NI'] = int(args[0])
#           Alternative values for the SIF file parameters:
# IE NI                  50             $-PARAMETER n=153, m=100
# IE NI                  100            $-PARAMETER n=303, m=200
# IE NI                  500            $-PARAMETER n=1503, m=1000
# IE NI                  1000           $-PARAMETER n=3003, m=2000 original value
# IE NI                  2000           $-PARAMETER n=6003, m=4000
# IE NI                  5000           $-PARAMETER n=15003, m=10000
        if nargin<2:
            v_['ALPHA'] = float(350.0);  #  SIF file default value
        else:
            v_['ALPHA'] = float(args[1])
        v_['RNI'] = float(v_['NI'])
        v_['NI-1'] = -1+v_['NI']
        v_['H'] = 1.0/v_['RNI']
        v_['H/4'] = 0.25*v_['H']
        v_['H/2'] = 0.5*v_['H']
        v_['AH'] = v_['ALPHA']*v_['H']
        v_['AH/2'] = 0.5*v_['AH']
        v_['-H/2'] = -0.5*v_['H']
        v_['0'] = 0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['0']),int(v_['NI'])+1):
            [iv,ix_,_] = s2x_ii('T'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'T'+str(I))
        for I in range(int(v_['0']),int(v_['NI'])+1):
            [iv,ix_,_] = s2x_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
        for I in range(int(v_['0']),int(v_['NI'])+1):
            [iv,ix_,_] = s2x_ii('U'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'U'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2x_ii('ENERGY',ig_)
        gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['0']),int(v_['NI-1'])+1):
            v_['I+1'] = 1+I
            [ig,ig_,_] = s2x_ii('EX'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'EX'+str(I))
            iv = ix_['X'+str(int(v_['I+1']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('ET'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'ET'+str(I))
            iv = ix_['T'+str(int(v_['I+1']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['T'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            iv = ix_['U'+str(int(v_['I+1']))]
            pbm.A[ig,iv] = float(v_['-H/2'])+pbm.A[ig,iv]
            iv = ix_['U'+str(I)]
            pbm.A[ig,iv] = float(v_['-H/2'])+pbm.A[ig,iv]
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
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        for I in range(int(v_['0']),int(v_['NI'])+1):
            pb.xlower[ix_['X'+str(I)]] = -0.05
            pb.xupper[ix_['X'+str(I)]] = 0.05
        for I in range(int(v_['0']),int(v_['NI'])+1):
            pb.xlower[ix_['T'+str(I)]] = -1.0
            pb.xupper[ix_['T'+str(I)]] = 1.0
        pb.xlower[ix_['X'+str(int(v_['0']))]] = 0.0
        pb.xupper[ix_['X'+str(int(v_['0']))]] = 0.0
        pb.xlower[ix_['X'+str(int(v_['NI']))]] = 0.0
        pb.xupper[ix_['X'+str(int(v_['NI']))]] = 0.0
        pb.xlower[ix_['T'+str(int(v_['0']))]] = 0.0
        pb.xupper[ix_['T'+str(int(v_['0']))]] = 0.0
        pb.xlower[ix_['T'+str(int(v_['NI']))]] = 0.0
        pb.xupper[ix_['T'+str(int(v_['NI']))]] = 0.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        for I in range(int(v_['0']),int(v_['NI'])+1):
            v_['RI'] = float(I)
            v_['TT'] = v_['RI']*v_['H']
            v_['CTT'] = np.cos(v_['TT'])
            v_['SCTT'] = 0.05*v_['CTT']
            pb.x0[ix_['X'+str(I)]] = float(v_['SCTT'])
            pb.x0[ix_['T'+str(I)]] = float(v_['SCTT'])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eCOS', iet_)
        elftv = loaset(elftv,it,0,'T')
        [it,iet_,_] = s2x_ii( 'eSIN', iet_)
        elftv = loaset(elftv,it,0,'T')
        [it,iet_,_] = s2x_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'U')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['0']),int(v_['NI'])+1):
            ename = 'C'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eCOS')
            ielftype = arrset(ielftype, ie, iet_["eCOS"])
            vname = 'T'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='T')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'S'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eSIN')
            ielftype = arrset(ielftype, ie, iet_["eSIN"])
            vname = 'T'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='T')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'USQ'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
            ielftype = arrset(ielftype, ie, iet_["eSQ"])
            vname = 'U'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='U')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['0']),int(v_['NI-1'])+1):
            v_['I+1'] = 1+I
            ig = ig_['EX'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['S'+str(int(v_['I+1']))])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-H/2']))
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['S'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-H/2']))
            ig = ig_['ENERGY']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['USQ'+str(int(v_['I+1']))])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['H/2']))
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['USQ'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['H/2']))
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C'+str(int(v_['I+1']))])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['AH/2']))
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['AH/2']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
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
        pb.pbclass = "OOR2-MN-V-V"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSQ(pbm,nargout,*args):

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

    @staticmethod
    def eCOS(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        CC = np.cos(EV_[0])
        SS = np.sin(EV_[0])
        f_   = CC
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -SS
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = -CC
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eSIN(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        CC = np.cos(EV_[0])
        SS = np.sin(EV_[0])
        f_   = SS
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = CC
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = -SS
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

