from s2mpjlib import *
class  DISCS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DISCS
#    *********
# 
#    The problem is to replace the nodes of a planar graph by discs
#    such that adjacent nodes have touching discs and non adjacent nodes
#    disjoint discs.  The smallest disc is constrained to have radius equal to
#    one and is centered at the origin.  The second discs is also
#    constrained to have its Y = 0.0, in order to cancel invariance under
#    rotation.
# 
#    Source:
#    W. Pulleyblank,
#    private communication, 1991.
# 
#    classification = "LQR2-MY-36-66"
# 
#    SIF input: A.R. Conn and Ph. Toint, April 1991.
# 
#    Number of nodes
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'DISCS'

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
        v_['NNODES'] = 12
        v_['EPSIL'] = 0.0001
        v_['1'] = 1
        v_['2'] = 2
        for I in range(int(v_['1']),int(v_['NNODES'])+1):
            v_['I-1'] = -1+I
            for J in range(int(v_['1']),int(v_['I-1'])+1):
                v_['A'+str(J)+','+str(I)] = 0.0
        v_['A1,2'] = 1.01
        v_['A1,7'] = 1.01
        v_['A2,3'] = 1.01
        v_['A2,4'] = 1.01
        v_['A3,4'] = 1.01
        v_['A4,5'] = 1.01
        v_['A4,6'] = 1.01
        v_['A5,6'] = 1.01
        v_['A5,11'] = 1.01
        v_['A6,7'] = 1.01
        v_['A7,8'] = 1.01
        v_['A7,9'] = 1.01
        v_['A8,9'] = 1.01
        v_['A8,10'] = 1.01
        v_['A9,10'] = 1.01
        v_['A10,11'] = 1.01
        v_['A10,12'] = 1.01
        v_['A11,12'] = 1.01
        v_['RNODES'] = float(v_['NNODES'])
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['NNODES'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
            [iv,ix_,_] = s2mpj_ii('Y'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'Y'+str(I))
            [iv,ix_,_] = s2mpj_ii('R'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'R'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['NNODES'])+1):
            [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['R'+str(I)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        for I in range(int(v_['2']),int(v_['NNODES'])+1):
            v_['I-1'] = -1+I
            for J in range(int(v_['1']),int(v_['I-1'])+1):
                v_['RAIJ'] = v_['A'+str(J)+','+str(I)]
                v_['AIJ'] = int(np.fix(v_['RAIJ']))
                for K in range(int(v_['1']),int(v_['AIJ'])+1):
                    [ig,ig_,_] = s2mpj_ii('S'+str(J)+','+str(I),ig_)
                    gtype = arrset(gtype,ig,'==')
                    cnames = arrset(cnames,ig,'S'+str(J)+','+str(I))
                v_['AIJ-1'] = -1+v_['AIJ']
                v_['NAIJ'] = -1*v_['AIJ-1']
                for K in range(int(v_['1']),int(v_['NAIJ'])+1):
                    [ig,ig_,_] = s2mpj_ii('S'+str(J)+','+str(I),ig_)
                    gtype = arrset(gtype,ig,'<=')
                    cnames = arrset(cnames,ig,'S'+str(J)+','+str(I))
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
        v_['-EPSIL'] = -1.0*v_['EPSIL']
        for I in range(int(v_['1']),int(v_['NNODES'])+1):
            v_['I-1'] = -1+I
            for J in range(int(v_['1']),int(v_['I-1'])+1):
                v_['RAIJ'] = v_['A'+str(J)+','+str(I)]
                v_['AIJ'] = int(np.fix(v_['RAIJ']))
                v_['AIJ-1'] = -1+v_['AIJ']
                v_['NAIJ'] = -1*v_['AIJ-1']
                for K in range(int(v_['1']),int(v_['NAIJ'])+1):
                    pbm.gconst  = (
                          arrset(pbm.gconst,ig_['S'+str(J)+','+str(I)],float(v_['-EPSIL'])))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower = np.zeros((pb.n,1))
        for I in range(int(v_['1']),int(v_['NNODES'])+1):
            pb.xlower[ix_['R'+str(I)]] = 1.0
        pb.xlower[ix_['X'+str(int(v_['1']))]] = 0.0
        pb.xupper[ix_['X'+str(int(v_['1']))]] = 0.0
        pb.xlower[ix_['Y'+str(int(v_['1']))]] = 0.0
        pb.xupper[ix_['Y'+str(int(v_['1']))]] = 0.0
        pb.xlower[ix_['Y'+str(int(v_['2']))]] = 0.0
        pb.xupper[ix_['Y'+str(int(v_['2']))]] = 0.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        for I in range(int(v_['1']),int(v_['NNODES'])+1):
            v_['RI'] = float(I)
            v_['RM'] = 0.03*v_['RI']
            v_['RP'] = 0.0055555*v_['RI']
            pb.x0[ix_['R'+str(I)]] = float(v_['RI'])
            pb.x0[ix_['X'+str(I)]] = float(v_['RM'])
            pb.x0[ix_['Y'+str(I)]] = float(v_['RP'])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eIPSQ', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        [it,iet_,_] = s2mpj_ii( 'eIMSQ', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['2']),int(v_['NNODES'])+1):
            v_['I-1'] = -1+I
            for J in range(int(v_['1']),int(v_['I-1'])+1):
                ename = 'DR'+str(J)+','+str(I)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eIPSQ')
                ielftype = arrset(ielftype, ie, iet_["eIPSQ"])
                vname = 'R'+str(I)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='X')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'R'+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                ename = 'DX'+str(J)+','+str(I)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eIMSQ')
                ielftype = arrset(ielftype, ie, iet_["eIMSQ"])
                vname = 'X'+str(I)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='X')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                ename = 'DY'+str(J)+','+str(I)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eIMSQ')
                ielftype = arrset(ielftype, ie, iet_["eIMSQ"])
                vname = 'Y'+str(I)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='X')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'Y'+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['2']),int(v_['NNODES'])+1):
            v_['I-1'] = -1+I
            for J in range(int(v_['1']),int(v_['I-1'])+1):
                ig = ig_['S'+str(J)+','+str(I)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['DR'+str(J)+','+str(I)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['DX'+str(J)+','+str(I)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
                posel = posel+1
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['DY'+str(J)+','+str(I)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
# ZL DISCS                              RNODES
#    Solution
# LO SOLTN(12)           20.46122911
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.cupper[np.arange(pb.nle)] = np.zeros((pb.nle,1))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "LQR2-MY-36-66"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eIPSQ(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((1,2))
        IV_ = np.zeros(1)
        U_[0,0] = U_[0,0]+1
        U_[0,1] = U_[0,1]+1
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
    def eIMSQ(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((1,2))
        IV_ = np.zeros(1)
        U_[0,0] = U_[0,0]+1
        U_[0,1] = U_[0,1]-1
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

