from s2mpjlib import *
class  DNIEPER(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DNIEPER
#    *********
# 
#    This problem models the planning of systematic use of water resources
#    in the basin of the river Dnieper.
# 
#    Source: p. 139sq in 
#    B.N. Pshenichnyj
#    "The Linearization Method for Constrained Optimization",
#    Springer Verlag, SCM Series 22, Heidelberg, 1994
# 
#    SIF input: Ph. Toint, December 1994.
# 
#    classification = "QOR2-MN-61-24"
# 
#    Problem data
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'DNIEPER'

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
        v_['C1'] = 5.61
        v_['C2'] = 4.68
        v_['C3'] = 1.62
        v_['C4'] = 1.8
        v_['C5'] = 2.13
        v_['C6'] = 2.1
        v_['C7'] = 1.99
        v_['C8'] = 2.02
        v_['C9'] = 2.14
        v_['C10'] = 2.15
        v_['C11'] = 2.36
        v_['C12'] = 2.63
        v_['C13'] = -0.02
        v_['C14'] = -0.01
        v_['C15'] = -0.16
        v_['C16'] = -0.47
        v_['C17'] = -0.75
        v_['C18'] = -0.94
        v_['C19'] = -0.93
        v_['C20'] = -0.99
        v_['C21'] = -0.42
        v_['C22'] = -0.07
        v_['C23'] = 0.04
        v_['C24'] = -0.06
        v_['1'] = 1
        v_['2'] = 2
        v_['4'] = 4
        v_['5'] = 5
        v_['8'] = 8
        v_['9'] = 9
        v_['12'] = 12
        v_['13'] = 13
        v_['14'] = 14
        v_['16'] = 16
        v_['17'] = 17
        v_['20'] = 20
        v_['21'] = 21
        v_['24'] = 24
        v_['25'] = 25
        v_['36'] = 36
        v_['37'] = 37
        v_['48'] = 48
        v_['49'] = 49
        v_['52'] = 52
        v_['53'] = 53
        v_['56'] = 56
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['56'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
        [iv,ix_,_] = s2mpj_ii('X0F',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X0F')
        [iv,ix_,_] = s2mpj_ii('X24F',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X24F')
        [iv,ix_,_] = s2mpj_ii('X12F',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X12F')
        [iv,ix_,_] = s2mpj_ii('X36F',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X36F')
        [iv,ix_,_] = s2mpj_ii('AC',ix_)
        pb.xnames=arrset(pb.xnames,iv,'AC')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['12'])+1):
            v_['I0'] = 12+I
            v_['I1'] = 24+I
            v_['I2'] = 36+I
            [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(int(v_['I1']))]
            pbm.A[ig,iv] = float(19.95)+pbm.A[ig,iv]
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = float(0.07656)+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['I2']))]
            pbm.A[ig,iv] = float(-24.89)+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['I0']))]
            pbm.A[ig,iv] = float(-0.7135)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        pbm.gscale = arrset(pbm.gscale,ig,float(-1.0))
        for I in range(int(v_['1']),int(v_['4'])+1):
            v_['I0'] = 24+I
            [ig,ig_,_] = s2mpj_ii('CC'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'CC'+str(I))
            iv = ix_['X'+str(int(v_['I0']))]
            pbm.A[ig,iv] = float(-2.68)+pbm.A[ig,iv]
        for I in range(int(v_['5']),int(v_['8'])+1):
            v_['I0'] = 24+I
            v_['I1'] = 44+I
            [ig,ig_,_] = s2mpj_ii('CC'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'CC'+str(I))
            iv = ix_['X'+str(int(v_['I0']))]
            pbm.A[ig,iv] = float(-2.68)+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['I1']))]
            pbm.A[ig,iv] = float(-2.68)+pbm.A[ig,iv]
        for I in range(int(v_['9']),int(v_['12'])+1):
            v_['I0'] = 24+I
            [ig,ig_,_] = s2mpj_ii('CC'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'CC'+str(I))
            iv = ix_['X'+str(int(v_['I0']))]
            pbm.A[ig,iv] = float(-2.68)+pbm.A[ig,iv]
        for I in range(int(v_['13']),int(v_['16'])+1):
            v_['I0'] = 12+I
            v_['I1'] = 24+I
            [ig,ig_,_] = s2mpj_ii('CC'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'CC'+str(I))
            iv = ix_['X'+str(int(v_['I0']))]
            pbm.A[ig,iv] = float(-2.68)+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['I1']))]
            pbm.A[ig,iv] = float(-2.68)+pbm.A[ig,iv]
            iv = ix_['AC']
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        for I in range(int(v_['17']),int(v_['20'])+1):
            v_['I0'] = 12+I
            v_['I1'] = 24+I
            v_['I2'] = 36+I
            [ig,ig_,_] = s2mpj_ii('CC'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'CC'+str(I))
            iv = ix_['X'+str(int(v_['I0']))]
            pbm.A[ig,iv] = float(-2.68)+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['I1']))]
            pbm.A[ig,iv] = float(-2.68)+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['I2']))]
            pbm.A[ig,iv] = float(-2.68)+pbm.A[ig,iv]
            iv = ix_['AC']
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        for I in range(int(v_['21']),int(v_['24'])+1):
            v_['I0'] = 12+I
            v_['I1'] = 24+I
            [ig,ig_,_] = s2mpj_ii('CC'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'CC'+str(I))
            iv = ix_['X'+str(int(v_['I0']))]
            pbm.A[ig,iv] = float(-2.68)+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['I1']))]
            pbm.A[ig,iv] = float(-2.68)+pbm.A[ig,iv]
            iv = ix_['AC']
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
        pbm.gconst = arrset(pbm.gconst,ig_['OBJ'],float(-112.464))
        for I in range(int(v_['1']),int(v_['24'])+1):
            v_['CST'] = -1.0*v_['C'+str(I)]
            pbm.gconst = arrset(pbm.gconst,ig_['CC'+str(I)],float(v_['CST']))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),float('inf'))
        for I in range(int(v_['1']),int(v_['12'])+1):
            pb.xlower[ix_['X'+str(I)]] = 51.2
            pb.xupper[ix_['X'+str(I)]] = 51.4
        for I in range(int(v_['13']),int(v_['24'])+1):
            pb.xlower[ix_['X'+str(I)]] = 15.0
            pb.xupper[ix_['X'+str(I)]] = 16.1
        for I in range(int(v_['25']),int(v_['36'])+1):
            pb.xlower[ix_['X'+str(I)]] = 0.4
            pb.xupper[ix_['X'+str(I)]] = 4.6
        for I in range(int(v_['37']),int(v_['48'])+1):
            pb.xlower[ix_['X'+str(I)]] = 0.5
            pb.xupper[ix_['X'+str(I)]] = 4.8
        for I in range(int(v_['49']),int(v_['56'])+1):
            pb.xlower[ix_['X'+str(I)]] = 0.0
            pb.xupper[ix_['X'+str(I)]] = 0.7
        pb.xlower[ix_['X0F']] = 50.82
        pb.xupper[ix_['X0F']] = 50.82
        pb.xlower[ix_['X24F']] = 2.0
        pb.xupper[ix_['X24F']] = 2.0
        pb.xlower[ix_['X12F']] = 15.5
        pb.xupper[ix_['X12F']] = 15.5
        pb.xlower[ix_['X36F']] = 2.3
        pb.xupper[ix_['X36F']] = 2.3
        pb.xlower[ix_['AC']] = -float('Inf')
        pb.xupper[ix_['AC']] = +float('Inf')
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        for I in range(int(v_['1']),int(v_['12'])+1):
            pb.x0[ix_['X'+str(I)]] = float(51.35)
        for I in range(int(v_['13']),int(v_['24'])+1):
            pb.x0[ix_['X'+str(I)]] = float(15.5)
        for I in range(int(v_['25']),int(v_['36'])+1):
            pb.x0[ix_['X'+str(I)]] = float(2.5)
        for I in range(int(v_['37']),int(v_['48'])+1):
            pb.x0[ix_['X'+str(I)]] = float(2.6)
        for I in range(int(v_['49']),int(v_['56'])+1):
            pb.x0[ix_['X'+str(I)]] = float(0.3)
        pb.x0[ix_['X0F']] = float(50.82)
        pb.x0[ix_['X24F']] = float(2.0)
        pb.x0[ix_['X12F']] = float(15.5)
        pb.x0[ix_['X36F']] = float(2.3)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eWJ', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        [it,iet_,_] = s2mpj_ii( 'eWK', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        [it,iet_,_] = s2mpj_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        [it,iet_,_] = s2mpj_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['1']),int(v_['12'])+1):
            v_['I1'] = 12+I
            v_['I2'] = 36+I
            ename = 'E'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset(ielftype, ie, iet_["en2PR"])
            vname = 'X'+str(int(v_['I1']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I2']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'ACSQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
        ielftype = arrset(ielftype, ie, iet_["eSQ"])
        vname = 'AC'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['1']),int(v_['12'])+1):
            v_['I0'] = 24+I
            ename = 'W1'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eWJ')
            ielftype = arrset(ielftype, ie, iet_["eWJ"])
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I0']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'W2'+str(int(v_['1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eWJ')
        ielftype = arrset(ielftype, ie, iet_["eWJ"])
        ename = 'W2'+str(int(v_['1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X0F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'W2'+str(int(v_['1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X24F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['2']),int(v_['12'])+1):
            v_['I0'] = 23+I
            v_['I1'] = -1+I
            ename = 'W2'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eWJ')
            ielftype = arrset(ielftype, ie, iet_["eWJ"])
            vname = 'X'+str(int(v_['I1']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I0']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['13']),int(v_['24'])+1):
            v_['I1'] = 24+I
            ename = 'W1'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eWK')
            ielftype = arrset(ielftype, ie, iet_["eWK"])
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I1']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'W2'+str(int(v_['13']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eWK')
        ielftype = arrset(ielftype, ie, iet_["eWK"])
        ename = 'W2'+str(int(v_['13']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X12F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'W2'+str(int(v_['13']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X36F'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['14']),int(v_['24'])+1):
            v_['I0'] = 23+I
            v_['I1'] = -1+I
            ename = 'W2'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eWK')
            ielftype = arrset(ielftype, ie, iet_["eWK"])
            vname = 'X'+str(int(v_['I1']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I0']))
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
        for I in range(int(v_['1']),int(v_['12'])+1):
            ig = ig_['OBJ']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(2.155))
        ig = ig_['OBJ']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['ACSQ'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-2000.0))
        for I in range(int(v_['1']),int(v_['24'])+1):
            ig = ig_['CC'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['W1'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['W2'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLUTION             1.87439D+04
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
        pb.pbclass = "QOR2-MN-61-24"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eWJ(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        A1 = 34.547
        A2 = -0.55878
        A3 = 8.05339
        A4 = -0.02252
        A5 = -0.29316
        A6 = -0.013521
        A7 = 0.00042
        A8 = 0.00267
        A9 = 0.000281
        A10 = 0.0000032
        f_   = (A1+A2*EV_[0]+A3*EV_[1]+A4*EV_[0]**2+A5*EV_[0]*EV_[1]+A6*EV_[1]**2+
             A7*EV_[0]**3+A8*EV_[0]**2*EV_[1]+A9*EV_[0]*EV_[1]**2+A10*EV_[1]**3)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = (A2+2.0*A4*EV_[0]+A5*EV_[1]+3.0*A7*EV_[0]**2+2.0*A8*EV_[0]*EV_[1]+
                 A9*EV_[1]**2)
            g_[1] = (A3+A5*EV_[0]+2.0*A6*EV_[1]+A8*EV_[0]**2+2.0*A9*EV_[0]*EV_[1]+
                 3.0*A10*EV_[1]**2)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 2.0*A4+6.0*A7*EV_[0]+2.0*A8*EV_[1]
                H_[0,1] = A5+2.0*A8*EV_[0]+2.0*A9*EV_[1]
                H_[1,0] = H_[0,1]
                H_[1,1] = 2.0*A6+2.0*A9*EV_[0]+6.0*A10*EV_[1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eWK(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        A1 = 20.923
        A2 = -4.22088
        A3 = 1.42061
        A4 = -0.41040
        A5 = -0.15082
        A7 = -0.00826
        A8 = 0.00404
        A9 = 0.000168
        A10 = -0.000038
        f_   = (A1+A2*EV_[0]+A3*EV_[1]+A4*EV_[0]**2+A5*EV_[0]*EV_[1]+A7*EV_[0]**3+
             A8*EV_[0]**2*EV_[1]+A9*EV_[0]*EV_[1]**2+A10*EV_[1]**3)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = (A2+2.0*A4*EV_[0]+A5*EV_[1]+3.0*A7*EV_[0]**2+2.0*A8*EV_[0]*EV_[1]+
                 A9*EV_[1]**2)
            g_[1] = A3+A5*EV_[0]+A8*EV_[0]**2+2.0*A9*EV_[0]*EV_[1]+3.0*A10*EV_[1]**2
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 2.0*A4+6.0*A7*EV_[0]+2.0*A8*EV_[1]
                H_[0,1] = A5+2.0*A8*EV_[0]+2.0*A9*EV_[1]
                H_[1,0] = H_[0,1]
                H_[1,1] = 2.0*A9*EV_[0]+6.0*A10*EV_[1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

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
    def eSQ(pbm,nargout,*args):

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
            g_[0] = EV_[0]+EV_[0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

