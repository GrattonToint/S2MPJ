from s2mpjlib import *
class  HIMMELBK(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HIMMELBK
#    *********
# 
#    A blending problem for multi-component mixtures, by Paviani.
#    It has a linear objective and linear and nonlinear constraints.
# 
#    Compared to the problem specified in Himmelblau, the inequality
#    constraints have been removed, because, as stated in this source,
#    they impose that
#    X(1)=X(2)=X(3)=X(7)=X(9)=X(9)=X(13)=X(14)=X(15)=X(19)=X(20)=X(21)=0
#    which is clearly contradictory with the given solution(s).  As
#    there does not seem to be a natural way to correct this statement
#    without knowing more about the original problem, the troublesome
#    constraints have been removed.
# 
#    Source: from problem 20 in
#    D.H. Himmelblau,
#    "Applied nonlinear programming",
#    McGraw-Hill, New-York, 1972.
# 
#    SIF input: Ph. Toint, March 1991.
# 
#    classification = "LOR2-MN-24-14"
# 
#    Problem data
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HIMMELBK'

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
        v_['F'] = 142.22471
        v_['B1'] = 44.094
        v_['B2'] = 58.12
        v_['B3'] = 58.12
        v_['B4'] = 137.4
        v_['B5'] = 120.9
        v_['B6'] = 170.9
        v_['B7'] = 62.501
        v_['B8'] = 84.94
        v_['B9'] = 133.425
        v_['B10'] = 82.507
        v_['B11'] = 46.07
        v_['B12'] = 60.097
        v_['B13'] = 44.094
        v_['B14'] = 58.12
        v_['B15'] = 58.12
        v_['B16'] = 137.4
        v_['B17'] = 120.9
        v_['B18'] = 170.9
        v_['B19'] = 62.501
        v_['B20'] = 84.94
        v_['B21'] = 133.425
        v_['B22'] = 82.507
        v_['B23'] = 46.07
        v_['B24'] = 60.097
        v_['C1'] = 123.7
        v_['C2'] = 31.7
        v_['C3'] = 45.7
        v_['C4'] = 14.7
        v_['C5'] = 84.7
        v_['C6'] = 27.7
        v_['C7'] = 49.7
        v_['C8'] = 7.1
        v_['C9'] = 2.1
        v_['C10'] = 17.7
        v_['C11'] = 0.85
        v_['C12'] = 0.64
        v_['D1'] = 123.7
        v_['D2'] = 31.7
        v_['D3'] = 45.7
        v_['D4'] = 14.7
        v_['D5'] = 84.7
        v_['D6'] = 27.7
        v_['D7'] = 49.7
        v_['D8'] = 7.1
        v_['D9'] = 2.1
        v_['D10'] = 17.7
        v_['D11'] = 0.85
        v_['D12'] = 0.64
        v_['1'] = 1
        v_['12'] = 12
        v_['13'] = 13
        v_['24'] = 24
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        pb.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for K in range(int(v_['1']),int(v_['24'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(K),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(K))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(0.0693)+pbm.A[ig,iv]
        iv = ix_['X2']
        pbm.A[ig,iv] = float(0.0577)+pbm.A[ig,iv]
        iv = ix_['X3']
        pbm.A[ig,iv] = float(0.05)+pbm.A[ig,iv]
        iv = ix_['X4']
        pbm.A[ig,iv] = float(0.2)+pbm.A[ig,iv]
        iv = ix_['X5']
        pbm.A[ig,iv] = float(0.26)+pbm.A[ig,iv]
        iv = ix_['X6']
        pbm.A[ig,iv] = float(0.55)+pbm.A[ig,iv]
        iv = ix_['X7']
        pbm.A[ig,iv] = float(0.06)+pbm.A[ig,iv]
        iv = ix_['X8']
        pbm.A[ig,iv] = float(0.1)+pbm.A[ig,iv]
        iv = ix_['X9']
        pbm.A[ig,iv] = float(0.12)+pbm.A[ig,iv]
        iv = ix_['X10']
        pbm.A[ig,iv] = float(0.18)+pbm.A[ig,iv]
        iv = ix_['X11']
        pbm.A[ig,iv] = float(0.1)+pbm.A[ig,iv]
        iv = ix_['X12']
        pbm.A[ig,iv] = float(0.09)+pbm.A[ig,iv]
        iv = ix_['X13']
        pbm.A[ig,iv] = float(0.0693)+pbm.A[ig,iv]
        iv = ix_['X14']
        pbm.A[ig,iv] = float(0.0577)+pbm.A[ig,iv]
        iv = ix_['X15']
        pbm.A[ig,iv] = float(0.05)+pbm.A[ig,iv]
        iv = ix_['X16']
        pbm.A[ig,iv] = float(0.2)+pbm.A[ig,iv]
        iv = ix_['X17']
        pbm.A[ig,iv] = float(0.26)+pbm.A[ig,iv]
        iv = ix_['X18']
        pbm.A[ig,iv] = float(0.55)+pbm.A[ig,iv]
        iv = ix_['X19']
        pbm.A[ig,iv] = float(0.06)+pbm.A[ig,iv]
        iv = ix_['X20']
        pbm.A[ig,iv] = float(0.1)+pbm.A[ig,iv]
        iv = ix_['X21']
        pbm.A[ig,iv] = float(0.12)+pbm.A[ig,iv]
        iv = ix_['X22']
        pbm.A[ig,iv] = float(0.18)+pbm.A[ig,iv]
        iv = ix_['X23']
        pbm.A[ig,iv] = float(0.1)+pbm.A[ig,iv]
        iv = ix_['X24']
        pbm.A[ig,iv] = float(0.09)+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['12'])+1):
            [ig,ig_,_] = s2mpj_ii('CA'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'CA'+str(I))
        for I in range(int(v_['1']),int(v_['24'])+1):
            [ig,ig_,_] = s2mpj_ii('CA13',ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'CA13')
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['12'])+1):
            v_['I+12'] = 12+I
            v_['1/DI'] = 1.0/v_['D'+str(I)]
            [ig,ig_,_] = s2mpj_ii('CA14',ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'CA14')
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = float(v_['1/DI'])+pbm.A[ig,iv]
            v_['F/BI+12'] = v_['F']/v_['B'+str(int(v_['I+12']))]
            iv = ix_['X'+str(int(v_['I+12']))]
            pbm.A[ig,iv] = float(v_['F/BI+12'])+pbm.A[ig,iv]
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
        pbm.gconst = arrset(pbm.gconst,ig_['CA13'],float(1.0))
        pbm.gconst = arrset(pbm.gconst,ig_['CA14'],float(1.671))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(0.04))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['1']),int(v_['12'])+1):
            v_['I+12'] = 12+I
            for J in range(int(v_['1']),int(v_['12'])+1):
                ename = 'E'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
                ielftype = arrset(ielftype, ie, iet_["en2PR"])
                vname = 'X'+str(int(v_['I+12']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.04)
                posev = find(elftv[ielftype[ie]],lambda x:x=='X')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.04)
                posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            for J in range(int(v_['13']),int(v_['24'])+1):
                ename = 'E'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
                ielftype = arrset(ielftype, ie, iet_["en2PR"])
                vname = 'X'+str(I)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.04)
                posev = find(elftv[ielftype[ie]],lambda x:x=='X')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.04)
                posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['12'])+1):
            v_['I+12'] = 12+I
            for J in range(int(v_['1']),int(v_['12'])+1):
                v_['BI/BJ'] = v_['B'+str(I)]/v_['B'+str(J)]
                v_['40BI/BJ'] = 40.0*v_['BI/BJ']
                ig = ig_['CA'+str(I)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['40BI/BJ']))
            for J in range(int(v_['13']),int(v_['24'])+1):
                v_['B+/BJ'] = v_['B'+str(int(v_['I+12']))]/v_['B'+str(J)]
                v_['CB+/BJ'] = v_['C'+str(I)]*v_['B+/BJ']
                v_['-CB+/BJ'] = -1.0*v_['CB+/BJ']
                ig = ig_['CA'+str(I)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-CB+/BJ']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN                0.0893344
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
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
        pb.pbclass = "LOR2-MN-24-14"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

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

