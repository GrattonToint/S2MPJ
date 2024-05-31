from s2xlib import *
class  BLOCKQP1(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : BLOCKQP1
#    *********
# 
#    A non-convex quadratic program with some structure.
# 
#    The objective function is of the form
#       sum (i=1,n) x_i y_i + sum(j=1,b) z_j^2
#    There are n equality constraints of the form
#         x_i + y_i + sum (j=1,b) z_j = b   
#    There is an inequality constraint of the form
#       sum(i=1,n) x_i + y_i + sum(j=1,b) z_j >= b + 1
#    Finally, there are simple bounds
#          1 <= x_i, y_i <= 1    (i=1,n)
#          0 <= z_j <= 2         (j=1,b)
# 
#    SIF input: Nick Gould, June 1994
# 
#    classification = "QLR2-AN-V-V"
# 
#    The number of equality constraints
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'BLOCKQP1'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'BLOCKQP1'
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
        if nargin<2:
            v_['B'] = int(5);  #  SIF file default value
        else:
            v_['B'] = int(args[1])
        v_['1'] = 1
        v_['RB'] = v_['B']
        v_['RB+1'] = 1+v_['RB']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2x_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
            [iv,ix_,_] = s2x_ii('Y'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'Y'+str(I))
        for J in range(int(v_['1']),int(v_['B'])+1):
            [iv,ix_,_] = s2x_ii('Z'+str(J),ix_)
            pb.xnames=arrset(pb.xnames,iv,'Z'+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2x_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['1']),int(v_['N'])+1):
            [ig,ig_,_] = s2x_ii(I,ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,I)
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = 1.0+pbm.A[ig,iv]
            iv = ix_['Y'+str(I)]
            pbm.A[ig,iv] = 1.0+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('E'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'E'+str(I))
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = 1.0+pbm.A[ig,iv]
            iv = ix_['Y'+str(I)]
            pbm.A[ig,iv] = 1.0+pbm.A[ig,iv]
        for J in range(int(v_['1']),int(v_['B'])+1):
            [ig,ig_,_] = s2x_ii('I',ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'I')
            iv = ix_['Z'+str(J)]
            pbm.A[ig,iv] = 1.0+pbm.A[ig,iv]
            for I in range(int(v_['1']),int(v_['N'])+1):
                [ig,ig_,_] = s2x_ii('E'+str(I),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'E'+str(I))
                iv = ix_['Z'+str(J)]
                pbm.A[ig,iv] = 1.0+pbm.A[ig,iv]
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
        for I in range(int(v_['1']),int(v_['N'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['E'+str(I)],v_['RB'])
        pbm.gconst = arrset(pbm.gconst,ig_['I'],v_['RB+1'])
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('inf'))
        pb.xupper = np.full((pb.n,1),float('inf'))
        for I in range(int(v_['1']),int(v_['N'])+1):
            pb.xlower[ix_['X'+str(I)]] = -1.0
            pb.xupper[ix_['X'+str(I)]] = 1.0
            pb.xlower[ix_['Y'+str(I)]] = -1.0
            pb.xupper[ix_['Y'+str(I)]] = 1.0
        for J in range(int(v_['1']),int(v_['B'])+1):
            pb.xlower[ix_['Z'+str(J)]] = 0.0
            pb.xupper[ix_['Z'+str(J)]] = 2.0
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),0.5)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'SQ', iet_)
        elftv = loaset(elftv,it,0,'Z')
        [it,iet_,_] = s2x_ii( 'PROD', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['1']),int(v_['N'])+1):
            ename = 'P'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'PROD')
            ielftype = arrset(ielftype, ie, iet_["PROD"])
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,0.5)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'Y'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,0.5)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['1']),int(v_['B'])+1):
            ename = 'S'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'SQ')
            ielftype = arrset(ielftype, ie, iet_["SQ"])
            vname = 'Z'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,0.5)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
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
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        for J in range(int(v_['1']),int(v_['B'])+1):
            ig = ig_['OBJ']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['S'+str(J)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
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
        pb.pbclass = "QLR2-AN-V-V"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def SQ(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = 0.5*EV_[0]*EV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 1.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def PROD(pbm,nargout,*args):

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

