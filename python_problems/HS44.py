from s2xlib import *
class  HS44(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS44
#    *********
# 
#    Source: problem 44 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: Ph.L. Toint, October 1990.
# 
#    classification = "QLR2-AN-4-6"
# 
#    some useful parameters, including N, the number of variables.
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS44'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'HS44'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 4
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
        iv = ix_['X1']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X2']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X3']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('CON1',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'CON1')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X2']
        pbm.A[ig,iv] = float(-2.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('CON2',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'CON2')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(-4.0)+pbm.A[ig,iv]
        iv = ix_['X2']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('CON3',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'CON3')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(-3.0)+pbm.A[ig,iv]
        iv = ix_['X2']
        pbm.A[ig,iv] = float(-4.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('CON4',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'CON4')
        iv = ix_['X3']
        pbm.A[ig,iv] = float(-2.0)+pbm.A[ig,iv]
        iv = ix_['X4']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('CON5',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'CON5')
        iv = ix_['X3']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X4']
        pbm.A[ig,iv] = float(-2.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('CON6',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'CON6')
        iv = ix_['X3']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X4']
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
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        pbm.gconst = arrset(pbm.gconst,ig_['CON1'],float(-8.0))
        pbm.gconst = arrset(pbm.gconst,ig_['CON2'],float(-12.0))
        pbm.gconst = arrset(pbm.gconst,ig_['CON3'],float(-12.0))
        pbm.gconst = arrset(pbm.gconst,ig_['CON4'],float(-8.0))
        pbm.gconst = arrset(pbm.gconst,ig_['CON5'],float(-8.0))
        pbm.gconst = arrset(pbm.gconst,ig_['CON6'],float(-5.0))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(0.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        ename = 'E1'
        [ie,ie_,newelt] = s2x_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'X1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E2'
        [ie,ie_,newelt] = s2x_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'X1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E3'
        [ie,ie_,newelt] = s2x_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'X2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E4'
        [ie,ie_,newelt] = s2x_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset( ielftype,ie,iet_['en2PR'])
        vname = 'X2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
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
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E1'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E2'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E4'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E3'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle+pb.neq,pb.m)] = np.zeros((pb.nge,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "QLR2-AN-4-6"
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

