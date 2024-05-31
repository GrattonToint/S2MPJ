from s2xlib import *
class  COSHFUN(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : COSHFUN
#    *********
# 
#    A nonlinear minmax problem.
# 
#    Source:
#    K. Jonasson and K. Madsen,
#    "Corrected sequential linear programming for sparse
#    minimax optimization", Technical report, Institute for Numerical
#    Analysis, Technical U. of Denmark.
# 
#    SIF input: Nick Gould, October 1992.
# 
#    classification = "LOR2-AN-V-V"
# 
#   the number of functions
# 
#           Alternative values for the SIF file parameters:
# IE M                   3              $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'COSHFUN'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'COSHFUN'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['M'] = int(8);  #  SIF file default value
        else:
            v_['M'] = int(args[0])
# IE M                   14             $-PARAMETER
# IE M                   20             $-PARAMETER     original value
# IE M                   200            $-PARAMETER
# IE M                   2000           $-PARAMETER
        v_['N'] = 3*v_['M']
        v_['N-3'] = -3+v_['N']
        v_['N-5'] = -5+v_['N']
        v_['0'] = 0
        v_['1'] = 1
        v_['3'] = 3
        v_['6'] = 6
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2x_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
        [iv,ix_,_] = s2x_ii('F',ix_)
        pb.xnames=arrset(pb.xnames,iv,'F')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2x_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['F']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['M'])+1):
            [ig,ig_,_] = s2x_ii('C'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'C'+str(I))
            iv = ix_['F']
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('C1',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'C1')
        iv = ix_['X3']
        pbm.A[ig,iv] = float(-2.0)+pbm.A[ig,iv]
        iv = ix_['X6']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        for I in range(int(v_['6']),int(v_['N-3'])+1,int(v_['3'])):
            v_['I-5'] = -5+I
            v_['I+3'] = 3+I
            v_['I/3'] = int(np.fix(I/v_['3']))
            [ig,ig_,_] = s2x_ii('C'+str(int(v_['I/3'])),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'C'+str(int(v_['I/3'])))
            iv = ix_['X'+str(int(v_['I-5']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = float(-2.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('C'+str(int(v_['I/3'])),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'C'+str(int(v_['I/3'])))
            iv = ix_['X'+str(int(v_['I+3']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('C'+str(int(v_['M'])),ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'C'+str(int(v_['M'])))
        iv = ix_['X'+str(int(v_['N-5']))]
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X'+str(int(v_['N']))]
        pbm.A[ig,iv] = float(-2.0)+pbm.A[ig,iv]
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
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eSQR', iet_)
        elftv = loaset(elftv,it,0,'X')
        [it,iet_,_] = s2x_ii( 'eCOSH', iet_)
        elftv = loaset(elftv,it,0,'X')
        [it,iet_,_] = s2x_ii( 'ePROD', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['3']),int(v_['N'])+1,int(v_['3'])):
            v_['I-1'] = -1+I
            v_['I-2'] = -2+I
            v_['I/3'] = int(np.fix(I/v_['3']))
            ename = 'SQR'+str(int(v_['I/3']))
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eSQR')
            ielftype = arrset(ielftype, ie, iet_["eSQR"])
            ename = 'SQR'+str(int(v_['I/3']))
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pb.x0 = np.zeros((pb.n,1))
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'COSH'+str(int(v_['I/3']))
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eCOSH')
            ielftype = arrset(ielftype, ie, iet_["eCOSH"])
            ename = 'COSH'+str(int(v_['I/3']))
            [ie,ie_,_] = s2x_ii(ename,ie_)
            vname = 'X'+str(int(v_['I-1']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'PROD'+str(int(v_['I/3']))
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'ePROD')
            ielftype = arrset(ielftype, ie, iet_["ePROD"])
            ename = 'PROD'+str(int(v_['I/3']))
            [ie,ie_,_] = s2x_ii(ename,ie_)
            vname = 'X'+str(int(v_['I-2']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'PROD'+str(int(v_['I/3']))
            [ie,ie_,_] = s2x_ii(ename,ie_)
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['M'])+1):
            ig = ig_['C'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['SQR'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['COSH'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['PROD'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.cupper[np.arange(pb.nle)] = np.zeros((pb.nle,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "LOR2-AN-V-V"
        pb.x0          = np.zeros((pb.n,1))
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
            g_[0] = 2.0e+0*EV_[0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0e+0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eCOSH(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        COSHX = np.cosh(EV_[0])
        f_   = COSHX
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = np.sinh(EV_[0])
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = COSHX
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ePROD(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = 2.0e+0*EV_[0]*EV_[0]*EV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 4.0e+0*EV_[0]*EV_[1]
            g_[1] = 2.0e+0*EV_[0]*EV_[0]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 4.0e+0*EV_[1]
                H_[0,1] = 4.0e+0*EV_[0]
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

