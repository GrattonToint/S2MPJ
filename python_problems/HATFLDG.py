from s2xlib import *
class  HATFLDG(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HATFLDG
#    *********
# 
#    A test problem from the OPTIMA user manual.
# 
#    Source:
#    "The OPTIMA user manual (issue No.8, p. 49)",
#    Numerical Optimization Centre, Hatfield Polytechnic (UK), 1989.
# 
#    SIF input: Ph. Toint, May 1990.
# 
#    classification = "NOR2-AY-25-25"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HATFLDG'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'HATFLDG'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 25
        v_['1'] = 1
        v_['2'] = 2
        v_['13'] = 13
        v_['N-1'] = -1+v_['N']
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
        for I in range(int(v_['1']),int(v_['N'])+1):
            [ig,ig_,_] = s2x_ii('G'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'G'+str(I))
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['13']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
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
        #%%%%%%%%%%%%%%%%%%  CONSTANTS %%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.full((ngrp,1),-1.0)
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(1.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        [it,iet_,_] = s2x_ii( 'en2PRI', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        elftv = loaset(elftv,it,2,'Z')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        ename = 'A1'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
        ielftype = arrset(ielftype, ie, iet_["en2PR"])
        vname = 'X1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['2']),int(v_['N-1'])+1):
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            ename = 'A'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'en2PRI')
            ielftype = arrset(ielftype, ie, iet_["en2PRI"])
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,1.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I-1']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,1.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I+1']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,1.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'A'+str(int(v_['N']))
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
        ielftype = arrset(ielftype, ie, iet_["en2PR"])
        ename = 'A'+str(int(v_['N']))
        [ie,ie_,_] = s2x_ii(ename,ie_)
        vname = 'X'+str(int(v_['N-1']))
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'A'+str(int(v_['N']))
        [ie,ie_,_] = s2x_ii(ename,ie_)
        vname = 'X'+str(int(v_['N']))
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        ig = ig_['G'+str(int(v_['1']))]
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['A'+str(int(v_['1']))])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        for I in range(int(v_['2']),int(v_['N'])+1):
            ig = ig_['G'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['A'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
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
        pb.pbclass = "NOR2-AY-25-25"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def en2PR(pbm,nargout,*args):

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
    def en2PRI(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((2,3))
        IV_ = np.zeros(2)
        U_[0,0] = U_[0,0]+1
        U_[1,1] = U_[1,1]+1
        U_[1,2] = U_[1,2]-1
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        f_   = IV_[0]*IV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = IV_[1]
            g_[1] = IV_[0]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = 1.0
                H_[1,0] = H_[0,1]
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

