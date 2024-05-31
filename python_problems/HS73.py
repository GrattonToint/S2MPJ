from s2xlib import *
class  HS73(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS73
#    *********
# 
#    A cattle feed problem
# 
#    Source: problem 73 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: B Debarsy, Apr 1990.
# 
#    classification = "LOR2-MN-4-3"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS73'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'HS73'
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
        pbm.A[ig,iv] = float(24.55)+pbm.A[ig,iv]
        iv = ix_['X2']
        pbm.A[ig,iv] = float(26.75)+pbm.A[ig,iv]
        iv = ix_['X3']
        pbm.A[ig,iv] = float(39.0)+pbm.A[ig,iv]
        iv = ix_['X4']
        pbm.A[ig,iv] = float(40.50)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('C1',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C1')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(2.3)+pbm.A[ig,iv]
        iv = ix_['X2']
        pbm.A[ig,iv] = float(5.6)+pbm.A[ig,iv]
        iv = ix_['X3']
        pbm.A[ig,iv] = float(11.1)+pbm.A[ig,iv]
        iv = ix_['X4']
        pbm.A[ig,iv] = float(1.3)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('C2',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C2')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(12.0)+pbm.A[ig,iv]
        iv = ix_['X2']
        pbm.A[ig,iv] = float(11.9)+pbm.A[ig,iv]
        iv = ix_['X3']
        pbm.A[ig,iv] = float(41.8)+pbm.A[ig,iv]
        iv = ix_['X4']
        pbm.A[ig,iv] = float(52.1)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('C3',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C3')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X4']
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
        pbm.gconst = arrset(pbm.gconst,ig_['C1'],float(5.0))
        pbm.gconst = arrset(pbm.gconst,ig_['C2'],float(21.0))
        pbm.gconst = arrset(pbm.gconst,ig_['C3'],float(1.0))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(1.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eSQRT', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        elftv = loaset(elftv,it,2,'U')
        elftv = loaset(elftv,it,3,'V')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        ename = 'E1'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSQRT')
        ielftype = arrset(ielftype, ie, iet_["eSQRT"])
        vname = 'X1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,1.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        ig = ig_['C2']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E1'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.645))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.clower[np.arange(pb.nle+pb.neq,pb.m)] = np.zeros((pb.nge,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "LOR2-MN-4-3"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSQRT(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        VALUE = 0.28*EV_[0]**2+0.19*EV_[1]**2+20.5*EV_[2]**2+0.62*EV_[3]**2
        FVALUE = np.sqrt(VALUE)
        FVAL = FVALUE*VALUE
        f_   = FVALUE
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = .28*EV_[0]/FVALUE
            g_[1] = .19*EV_[1]/FVALUE
            g_[2] = 20.5*EV_[2]/FVALUE
            g_[3] = .62*EV_[3]/FVALUE
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,0] = -.0784*EV_[0]**2/FVAL+.28/FVALUE
                H_[0,1] = -.0532*EV_[0]*EV_[1]/FVAL
                H_[1,0] = H_[0,1]
                H_[0,2] = -5.74*EV_[0]*EV_[2]/FVAL
                H_[2,0] = H_[0,2]
                H_[0,3] = -.1736*EV_[0]*EV_[3]/FVAL
                H_[3,0] = H_[0,3]
                H_[1,1] = -.0361*EV_[1]**2/FVAL+.19/FVALUE
                H_[1,2] = -3.895*EV_[1]*EV_[2]/FVAL
                H_[2,1] = H_[1,2]
                H_[1,3] = -.1178*EV_[1]*EV_[3]/FVAL
                H_[3,1] = H_[1,3]
                H_[2,2] = -420.25*EV_[2]**2/FVAL+20.5/FVALUE
                H_[2,3] = -12.71*EV_[2]*EV_[3]/FVAL
                H_[3,2] = H_[2,3]
                H_[3,3] = -.3844*EV_[3]**2/FVAL+.62/FVALUE
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

