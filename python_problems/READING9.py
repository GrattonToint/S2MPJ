from s2mpjlib import *
class  READING9(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : READING9
#    *********
# 
#    A nonlinear optimal control problem from Nancy Nichols
#    with a given initial condition.
#    This problem arises in tide modelling.
# 
#    Source: a variant upon a problem in
#    S. Lyle and N.K. Nichols,
#    "Numerical Methods for Optimal Control Problems with State Constraints",
#    Numerical Analysis Report 8/91, Dept of Mathematics, 
#    University of Reading, UK.
# 
#    SIF input: Nick Gould and Ph. Toint, March 1995
# 
#    classification = "OOR2-MN-V-V"
# 
#    Number of discretized points in [0,1] - 1
# 
#           Alternative values for the SIF file parameters:
# IE N+1                 3              $-PARAMETER n=6, m=2
# IE N+1                 51             $-PARAMETER n=102, m=50
# IE N+1                 101            $-PARAMETER n=202, m=100
# IE N+1                 201            $-PARAMETER n=402, m=200
# IE N+1                 501            $-PARAMETER n=1002, m=500  original value
# IE N+1                 1001           $-PARAMETER n=2002, m=1000
# IE N+1                 5001           $-PARAMETER n=10002, m= 5000
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'READING9'

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
        if nargin<1:
            v_['N+1'] = int(3);  #  SIF file default value
        else:
            v_['N+1'] = int(args[0])
        v_['N'] = -1+v_['N+1']
        v_['N-1'] = -1+v_['N']
        v_['RN'] = float(v_['N'])
        v_['H'] = 1.0/v_['RN']
        v_['-H/2'] = -0.5*v_['H']
        v_['K1'] = 0.07716
        v_['K2'] = 0.07716
        v_['K1H'] = v_['K1']*v_['H']
        v_['K1H+1'] = 1.0+v_['K1H']
        v_['-K1H-1'] = -1.0*v_['K1H+1']
        v_['K2H'] = v_['K2']*v_['H']
        v_['0'] = 0
        v_['1'] = 1
        v_['2'] = 2
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['0']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('P'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'P'+str(I))
        for I in range(int(v_['0']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('U'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'U'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['0']),int(v_['N-1'])+1):
            v_['I+1'] = 1+I
            [ig,ig_,_] = s2mpj_ii('S'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'S'+str(I))
            iv = ix_['P'+str(int(v_['I+1']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['P'+str(I)]
            pbm.A[ig,iv] = float(v_['-K1H-1'])+pbm.A[ig,iv]
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
        for I in range(int(v_['0']),int(v_['N-1'])+1):
            v_['RI'] = float(I)
            v_['T'] = v_['RI']*v_['H']
            v_['SINT'] = np.sin(v_['T'])
            v_['HSINT'] = v_['H']*v_['SINT']
            pbm.gconst = arrset(pbm.gconst,ig_['S'+str(I)],float(v_['HSINT']))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),float('inf'))
        for I in range(int(v_['0']),int(v_['N'])+1):
            pb.xupper[ix_['U'+str(I)]] = 1.0
            pb.xlower[ix_['P'+str(I)]] = -float('Inf')
            pb.xupper[ix_['P'+str(I)]] = +float('Inf')
        pb.xlower[ix_['P'+str(int(v_['0']))]] = 0.0
        pb.xupper[ix_['P'+str(int(v_['0']))]] = 0.0
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(0.2))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'ePROD', iet_)
        elftv = loaset(elftv,it,0,'U')
        elftv = loaset(elftv,it,1,'P')
        [it,iet_,_] = s2mpj_ii( 'ePROD2', iet_)
        elftv = loaset(elftv,it,0,'U')
        elftv = loaset(elftv,it,1,'P')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['0']),int(v_['N'])+1):
            ename = 'OE'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'ePROD2')
            ielftype = arrset(ielftype, ie, iet_["ePROD2"])
            vname = 'P'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.2)
            posev = find(elftv[ielftype[ie]],lambda x:x=='P')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'U'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.2)
            posev = find(elftv[ielftype[ie]],lambda x:x=='U')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['0']),int(v_['N-1'])+1):
            ename = 'CE'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'ePROD')
            ielftype = arrset(ielftype, ie, iet_["ePROD"])
            vname = 'P'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.2)
            posev = find(elftv[ielftype[ie]],lambda x:x=='P')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'U'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.2)
            posev = find(elftv[ielftype[ie]],lambda x:x=='U')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['0']),int(v_['N-1'])+1):
            v_['I+1'] = 1+I
            ig = ig_['OBJ']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['OE'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-H/2']))
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['OE'+str(int(v_['I+1']))])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-H/2']))
            ig = ig_['S'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['CE'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['K2H']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# XL SOLUTION            -4.41677D-02   $ (n=500)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
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
# ********************
#  SET UP THE GROUPS *
#  ROUTINE           *
# ********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def ePROD(pbm,nargout,*args):

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

    @staticmethod
    def ePROD2(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = (EV_[0]*EV_[1])**2
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*EV_[0]*EV_[1]**2
            g_[1] = 2.0*EV_[1]*EV_[0]**2
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 2.0*EV_[1]**2
                H_[0,1] = 4.0*EV_[0]*EV_[1]
                H_[1,0] = H_[0,1]
                H_[1,1] = 2.0*EV_[0]**2
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

