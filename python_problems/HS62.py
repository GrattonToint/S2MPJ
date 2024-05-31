from s2xlib import *
class  HS62(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS62
#    *********
# 
#    Source: problem 62 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: J-M Collin and Ph. Toint, April 1990.
# 
#    classification = "OLR2-AY-3-1"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS62'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'HS62'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 3
        v_['1'] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2x_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2x_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2x_ii('C1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C1')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
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
        pbm.gconst = arrset(pbm.gconst,ig_['C1'],float(1.0))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('inf'))
        pb.xupper = np.full((pb.n,1),float('inf'))
        for I in range(int(v_['1']),int(v_['N'])+1):
            pb.xlower[ix_['X'+str(I)]] = 0.0
            pb.xupper[ix_['X'+str(I)]] = 1.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        if('X1' in ix_):
            pb.x0[ix_['X1']] = float(0.7)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X1']),float(0.7)))
        if('X2' in ix_):
            pb.x0[ix_['X2']] = float(0.2)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X2']),float(0.2)))
        if('X3' in ix_):
            pb.x0[ix_['X3']] = float(0.1)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X3']),float(0.1)))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eLN3N', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        elftv = loaset(elftv,it,2,'Z')
        [it,iet_,_] = s2x_ii( 'eLN3D', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        elftv = loaset(elftv,it,2,'Z')
        [it,iet_,_] = s2x_ii( 'eLN2N', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        [it,iet_,_] = s2x_ii( 'eLN2D', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        [it,iet_,_] = s2x_ii( 'eLN1N', iet_)
        elftv = loaset(elftv,it,0,'X')
        [it,iet_,_] = s2x_ii( 'eLN1D', iet_)
        elftv = loaset(elftv,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        ename = 'OE1N'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eLN3N')
        ielftype = arrset(ielftype, ie, iet_["eLN3N"])
        vname = 'X1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'OE1D'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eLN3D')
        ielftype = arrset(ielftype, ie, iet_["eLN3D"])
        vname = 'X1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'OE2N'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eLN2N')
        ielftype = arrset(ielftype, ie, iet_["eLN2N"])
        vname = 'X2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'OE2D'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eLN2D')
        ielftype = arrset(ielftype, ie, iet_["eLN2D"])
        vname = 'X2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'OE3N'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eLN1N')
        ielftype = arrset(ielftype, ie, iet_["eLN1N"])
        vname = 'X3'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'OE3D'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eLN1D')
        ielftype = arrset(ielftype, ie, iet_["eLN1D"])
        vname = 'X3'
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
        ig = ig_['OBJ']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['OE1N'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-8204.37))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['OE1D'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(8204.37))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['OE2N'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-9008.72))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['OE2D'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(9008.72))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['OE3N'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-9330.46))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['OE3D'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(9330.46))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
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
        pb.pbclass = "OLR2-AY-3-1"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eLN3N(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((1,3))
        IV_ = np.zeros(1)
        U_[0,0] = U_[0,0]+1
        U_[0,1] = U_[0,1]+1
        U_[0,2] = U_[0,2]+1
        IV_[0] = U_[0:1,:].dot(EV_)
        NUM = IV_[0]+0.03
        f_   = np.log(NUM)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 1.0/NUM
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = -1.0/NUM**2
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eLN3D(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((1,3))
        IV_ = np.zeros(1)
        U_[0,0] = U_[0,0]+9.000000e-02
        U_[0,1] = U_[0,1]+1
        U_[0,2] = U_[0,2]+1
        IV_[0] = U_[0:1,:].dot(EV_)
        NUM = IV_[0]+0.03
        f_   = np.log(NUM)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 1.0/NUM
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = -1.0/NUM**2
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eLN2N(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((1,2))
        IV_ = np.zeros(1)
        U_[0,0] = U_[0,0]+1
        U_[0,1] = U_[0,1]+1
        IV_[0] = U_[0:1,:].dot(EV_)
        NUM = IV_[0]+0.03
        f_   = np.log(NUM)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 1.0/NUM
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = -1.0/NUM**2
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eLN2D(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((1,2))
        IV_ = np.zeros(1)
        U_[0,0] = U_[0,0]+7.000000e-02
        U_[0,1] = U_[0,1]+1
        IV_[0] = U_[0:1,:].dot(EV_)
        NUM = IV_[0]+0.03
        f_   = np.log(NUM)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 1.0/NUM
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = -1.0/NUM**2
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eLN1N(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        NUM = EV_[0]+0.03
        f_   = np.log(NUM)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 1.0/NUM
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = -1.0/NUM**2
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eLN1D(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        NUM = 0.13*EV_[0]+0.03
        f_   = np.log(NUM)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 0.13/NUM
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = -0.0169/NUM**2
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

