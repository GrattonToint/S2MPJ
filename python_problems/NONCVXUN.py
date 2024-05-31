from s2xlib import *
class  NONCVXUN(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : NONCVXUN
#    *********
# 
#    A nonconvex unconstrained function with a unique minimum value
# 
#    SIF input: Nick Gould, April 1996
# 
#    classification = "OUR2-AN-V-0"
# 
#    The number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'NONCVXUN'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'NONCVXUN'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(10);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
#           Alternative values for the SIF file parameters:
# IE N                   100            $-PARAMETER
# IE N                   1000           $-PARAMETER     original value
# IE N                   5000           $-PARAMETER
# IE N                   10000          $-PARAMETER
# IE N                   100000         $-PARAMETER
        v_['1'] = 1
        v_['2'] = 2
        v_['4'] = 4
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
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['RI'] = float(I)
            pb.x0[ix_['X'+str(I)]] = float(v_['RI'])
        pass
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        [it,iet_,_] = s2x_ii( 'eCOS', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['1']),int(v_['N'])+1):
            ename = 'SQ'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
            ielftype = arrset(ielftype, ie, iet_["eSQ"])
            ename = 'COS'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eCOS')
            ielftype = arrset(ielftype, ie, iet_["eCOS"])
            ename = 'SQ'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'COS'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            v_['J'] = 2*I
            v_['J'] = -1+v_['J']
            v_['K'] = int(np.fix(v_['J']/v_['N']))
            v_['K'] = v_['K']*v_['N']
            v_['J'] = v_['J']-v_['K']
            v_['J'] = 1+v_['J']
            ename = 'SQ'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            vname = 'X'+str(int(v_['J']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'COS'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            vname = 'X'+str(int(v_['J']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            v_['J'] = 3*I
            v_['J'] = -1+v_['J']
            v_['K'] = int(np.fix(v_['J']/v_['N']))
            v_['K'] = v_['K']*v_['N']
            v_['J'] = v_['J']-v_['K']
            v_['J'] = 1+v_['J']
            ename = 'SQ'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            vname = 'X'+str(int(v_['J']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'COS'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            vname = 'X'+str(int(v_['J']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            ig = ig_['OBJ']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['SQ'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['COS'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "OUR2-AN-V-0"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSQ(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((1,3))
        IV_ = np.zeros(1)
        U_[0,0] = U_[0,0]+1
        U_[0,1] = U_[0,1]+1
        U_[0,2] = U_[0,2]+1
        IV_[0] = U_[0:1,:].dot(EV_)
        f_   = IV_[0]*IV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = IV_[0]+IV_[0]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0
                H_ = U_.T.dot(H_).dot(U_)
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
        U_ = np.zeros((1,3))
        IV_ = np.zeros(1)
        U_[0,0] = U_[0,0]+1
        U_[0,1] = U_[0,1]+1
        U_[0,2] = U_[0,2]+1
        IV_[0] = U_[0:1,:].dot(EV_)
        COSV = 4.0*np.cos(IV_[0])
        f_   = COSV
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -4.0*np.sin(IV_[0])
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = -COSV
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

