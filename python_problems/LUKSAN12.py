from s2xlib import *
class  LUKSAN12(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : LUKSAN12
#    *********
# 
#    Problem 12 (chained and modified HS47) in the paper
# 
#      L. Luksan
#      Hybrid methods in large sparse nonlinear least squares
#      J. Optimization Theory & Applications 89(3) 575-595 (1996)
# 
#    SIF input: Nick Gould, June 2017.
# 
#    classification = "NOR2-AN-V-V"
# 
#   seed for dimensions
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'LUKSAN12'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'LUKSAN12'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['S'] = 32
        v_['N'] = 3*v_['S']
        v_['N'] = 2+v_['N']
        v_['M'] = 6*v_['S']
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
        v_['I'] = 1
        v_['K'] = 1
        for J in range(int(v_['1']),int(v_['S'])+1):
            v_['K+1'] = 1+v_['K']
            v_['K+2'] = 2+v_['K']
            v_['K+3'] = 3+v_['K']
            v_['K+4'] = 4+v_['K']
            v_['K+5'] = 5+v_['K']
            v_['I+1'] = 1+v_['I']
            v_['I+2'] = 2+v_['I']
            [ig,ig_,_] = s2x_ii('E'+str(int(v_['K'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'E'+str(int(v_['K'])))
            iv = ix_['X'+str(int(v_['I+1']))]
            pbm.A[ig,iv] = float(-10.0e0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('E'+str(int(v_['K+1'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'E'+str(int(v_['K+1'])))
            iv = ix_['X'+str(int(v_['I+2']))]
            pbm.A[ig,iv] = float(1.0e0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('E'+str(int(v_['K+2'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'E'+str(int(v_['K+2'])))
            [ig,ig_,_] = s2x_ii('E'+str(int(v_['K+3'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'E'+str(int(v_['K+3'])))
            [ig,ig_,_] = s2x_ii('E'+str(int(v_['K+4'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'E'+str(int(v_['K+4'])))
            [ig,ig_,_] = s2x_ii('E'+str(int(v_['K+5'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'E'+str(int(v_['K+5'])))
            iv = ix_['X'+str(int(v_['I+1']))]
            pbm.A[ig,iv] = float(1.0e0)+pbm.A[ig,iv]
            v_['I'] = 3+v_['I']
            v_['K'] = 6+v_['K']
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
        v_['K'] = 1
        for J in range(int(v_['1']),int(v_['S'])+1):
            v_['K+1'] = 1+v_['K']
            v_['K+4'] = 4+v_['K']
            v_['K+5'] = 5+v_['K']
            pbm.gconst = arrset(pbm.gconst,ig_['E'+str(int(v_['K+1']))],float(1.0))
            pbm.gconst = arrset(pbm.gconst,ig_['E'+str(int(v_['K+4']))],float(10.0))
            pbm.gconst = arrset(pbm.gconst,ig_['E'+str(int(v_['K+5']))],float(20.0))
            v_['K'] = 6+v_['K']
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        for I in range(int(v_['1']),int(v_['N'])+1):
            pb.x0[ix_['X'+str(I)]] = float(-1.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eE1', iet_)
        elftv = loaset(elftv,it,0,'X0')
        [it,iet_,_] = s2x_ii( 'eE3', iet_)
        elftv = loaset(elftv,it,0,'X3')
        [it,iet_,_] = s2x_ii( 'eE4', iet_)
        elftv = loaset(elftv,it,0,'X4')
        [it,iet_,_] = s2x_ii( 'eE5', iet_)
        elftv = loaset(elftv,it,0,'X0')
        elftv = loaset(elftv,it,1,'X3')
        [it,iet_,_] = s2x_ii( 'eF5', iet_)
        elftv = loaset(elftv,it,0,'X3')
        elftv = loaset(elftv,it,1,'X4')
        [it,iet_,_] = s2x_ii( 'eE6', iet_)
        elftv = loaset(elftv,it,0,'X2')
        elftv = loaset(elftv,it,1,'X3')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        v_['I'] = 1
        v_['K'] = 1
        for J in range(int(v_['1']),int(v_['S'])+1):
            v_['K+2'] = 2+v_['K']
            v_['K+3'] = 3+v_['K']
            v_['K+4'] = 4+v_['K']
            v_['K+5'] = 5+v_['K']
            v_['I+2'] = 2+v_['I']
            v_['I+3'] = 3+v_['I']
            v_['I+4'] = 4+v_['I']
            ename = 'E'+str(int(v_['K']))
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eE1')
            ielftype = arrset(ielftype, ie, iet_["eE1"])
            ename = 'E'+str(int(v_['K']))
            [ie,ie_,_] = s2x_ii(ename,ie_)
            vname = 'X'+str(int(v_['I']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X0')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'E'+str(int(v_['K+2']))
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eE3')
            ielftype = arrset(ielftype, ie, iet_["eE3"])
            ename = 'E'+str(int(v_['K+2']))
            [ie,ie_,_] = s2x_ii(ename,ie_)
            vname = 'X'+str(int(v_['I+3']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'E'+str(int(v_['K+3']))
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eE4')
            ielftype = arrset(ielftype, ie, iet_["eE4"])
            ename = 'E'+str(int(v_['K+3']))
            [ie,ie_,_] = s2x_ii(ename,ie_)
            vname = 'X'+str(int(v_['I+4']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X4')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'E'+str(int(v_['K+4']))
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eE5')
            ielftype = arrset(ielftype, ie, iet_["eE5"])
            ename = 'E'+str(int(v_['K+4']))
            [ie,ie_,_] = s2x_ii(ename,ie_)
            vname = 'X'+str(int(v_['I']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X0')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'E'+str(int(v_['K+4']))
            [ie,ie_,_] = s2x_ii(ename,ie_)
            vname = 'X'+str(int(v_['I+3']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'F'+str(int(v_['K+4']))
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eF5')
            ielftype = arrset(ielftype, ie, iet_["eF5"])
            ename = 'F'+str(int(v_['K+4']))
            [ie,ie_,_] = s2x_ii(ename,ie_)
            vname = 'X'+str(int(v_['I+3']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'F'+str(int(v_['K+4']))
            [ie,ie_,_] = s2x_ii(ename,ie_)
            vname = 'X'+str(int(v_['I+4']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X4')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'E'+str(int(v_['K+5']))
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eE6')
            ielftype = arrset(ielftype, ie, iet_["eE6"])
            ename = 'E'+str(int(v_['K+5']))
            [ie,ie_,_] = s2x_ii(ename,ie_)
            vname = 'X'+str(int(v_['I+2']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'E'+str(int(v_['K+5']))
            [ie,ie_,_] = s2x_ii(ename,ie_)
            vname = 'X'+str(int(v_['I+3']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            v_['I'] = 3+v_['I']
            v_['K'] = 6+v_['K']
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        v_['K'] = 1
        for J in range(int(v_['1']),int(v_['S'])+1):
            v_['K+2'] = 2+v_['K']
            v_['K+3'] = 3+v_['K']
            v_['K+4'] = 4+v_['K']
            v_['K+5'] = 5+v_['K']
            ig = ig_['E'+str(int(v_['K']))]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['K']))])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            ig = ig_['E'+str(int(v_['K+2']))]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['K+2']))])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            ig = ig_['E'+str(int(v_['K+3']))]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['K+3']))])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            ig = ig_['E'+str(int(v_['K+4']))]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['K+4']))])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F'+str(int(v_['K+4']))])
            pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
            ig = ig_['E'+str(int(v_['K+5']))]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['K+5']))])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            v_['K'] = 6+v_['K']
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
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
        pb.pbclass = "NOR2-AN-V-V"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eE1(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = 10.0e0*EV_[0]**2
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 20.0e0*EV_[0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 20.0e0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eE3(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = (EV_[0]-1.0e0)**2
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0e0*(EV_[0]-1.0e0)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0e0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eE4(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = (EV_[0]-1.0e0)**3
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 3.0e0*(EV_[0]-1.0e0)**2
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 6.0e0*(EV_[0]-1.0e0)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eE5(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[1]*EV_[0]*EV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0e0*EV_[1]*EV_[0]
            g_[1] = EV_[0]*EV_[0]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 2.0e0*EV_[1]
                H_[0,1] = 2.0e0*EV_[0]
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eF5(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((1,2))
        IV_ = np.zeros(1)
        U_[0,0] = U_[0,0]+1
        U_[0,1] = U_[0,1]-1
        IV_[0] = U_[0:1,:].dot(EV_)
        f_   = np.sin(IV_[0])
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = np.cos(IV_[0])
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = -np.sin(IV_[0])
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eE6(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = (EV_[0]**4)*(EV_[1]**2)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 4.0e0*(EV_[0]**3)*(EV_[1]**2)
            g_[1] = 2.0e0*(EV_[0]**4)*EV_[1]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 12.0e0*(EV_[0]**2)*(EV_[1]**2)
                H_[0,1] = 8.0e0*(EV_[0]**3)*EV_[1]
                H_[1,0] = H_[0,1]
                H_[1,1] = 2.0e0*(EV_[0]**4)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
